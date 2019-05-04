package main

//
// Merge 2 to n VCF files covering the same genomic range
//-------------------------------------------------------
// The are three major aspects:
// 1) Handling file merge as an extension of a two-way merge - this can be regarded
// as an implemetation of 'direct k-way' merge operating on files
//
// 2) For minimum-key records sets record combination takes place, if set-size > 1
//
// args:
//  --tpltfile: a text file of template file paths for files to be merged
//  --paramfile: file of parameters for genotype resolution
//  --chr: chromosome
//  --logfile: full filepath for logging
//  --vcfprfx: directory root for vcf files
//
import (
	"bufio"
	"compress/gzip"
	"flag"
	"fmt"
	"genometrics"
	"io"
	"log"
	"os"
	"sample"
	"sort"
	"strings"
	"variant"
	"vcfmerge"
)

const max_posn int64 = 9999999999

var empty_record = []string{}

//-----------------------------------------------
// global vars, accessed by multiple funcs
//-----------------------------------------------
var tpltFilePath string
var paramFilePath string
var logFilePath string
var vcfPathPref string
var chr string
var threshold float64

//-----------------------------------------------
// main package routines
//-----------------------------------------------
func init() {
	const (
		defaultTpltFilePath  = "./data/vcf_file_template.txt"
		tusage               = "File template strings"
		defaultParamFilePath = "./data/params.cfg"
		pusage               = "QC Parameter file"
		defaultLogFilePath   = "./data/filemergevcf_output.log"
		lusage               = "Log file"
		defaultvcfPathPref   = "/var/data"
		vusage               = "default path prefix for vcf files"
		defaultThreshold     = 0.9
		thrusage             = "Prob threshold"
		defaultChr           = "22"
		chrusage             = "default chromosome (number as string)"
	)
	flag.StringVar(&tpltFilePath, "tpltfile", defaultTpltFilePath, tusage)
	flag.StringVar(&tpltFilePath, "t", defaultTpltFilePath, tusage+" (shorthand)")
	flag.StringVar(&paramFilePath, "paramfile", defaultParamFilePath, pusage)
	flag.StringVar(&paramFilePath, "p", defaultParamFilePath, pusage+" (shorthand)")
	flag.StringVar(&logFilePath, "logfile", defaultLogFilePath, lusage)
	flag.StringVar(&logFilePath, "l", defaultLogFilePath, lusage+" (shorthand)")
	flag.StringVar(&vcfPathPref, "vcfprfx", defaultvcfPathPref, vusage)
	flag.StringVar(&vcfPathPref, "v", defaultvcfPathPref, vusage+" (shorthand)")
	flag.Float64Var(&threshold, "threshold", defaultThreshold, thrusage)
	flag.Float64Var(&threshold, "h", defaultThreshold, thrusage+" (shorthand)")
	flag.StringVar(&chr, "chr", defaultChr, chrusage)
	flag.StringVar(&chr, "c", defaultChr, chrusage+" (shorthand)")
	flag.Parse()
}

func check(e error) {
	if e != nil {
		log.Println("err != nil")
		log.Fatal(e)
	}
}

func main() {
	// set up logging to a file
	lf, err := os.OpenFile(logFilePath, os.O_RDWR|os.O_CREATE|os.O_APPEND, 0666)
	check(err)
	defer lf.Close()

	log.SetOutput(lf)
	log.Printf("START merge %s, %s\n", paramFilePath, tpltFilePath)

	// Load file templates
	f, err := os.Open(tpltFilePath)
	check(err)
	defer f.Close()
	scanner := bufio.NewScanner(f)
	assaytype_filename := make(map[string]string)
	assaytype_list := make([]string, 0)
	fscanners := make(map[string]*bufio.Scanner)
	freaders := make(map[string]*bufio.Reader)

	for scanner.Scan() {
		text := scanner.Text()
		if !strings.HasPrefix(text, "#") {
			fields := strings.Split(text, "=")
			assaytype_filename[fields[0]] = fmt.Sprintf(fields[1], vcfPathPref, chr)
			assaytype_list = append(assaytype_list, fields[0])
		}
	}
	// Load QC parameters (if present)
	fp, err := os.Open(paramFilePath)
	//check(err)
	defer fp.Close()
	testnum := ""
	callrate := ""
	mafdelta := ""
	infoscore := ""
	scanner = bufio.NewScanner(fp)
	for scanner.Scan() {
		text := scanner.Text()
		if !strings.HasPrefix(text, "#") {
			fields := strings.Split(text, "=")
			if fields[0] == "TESTNUM" {
				testnum = fields[1]
			}
			if fields[0] == "CALLRATE" {
				callrate = fields[1]
			}
			if fields[0] == "MAFDELTA" {
				mafdelta = fields[1]
			}
			if fields[0] == "INFOSCORE" {
				infoscore = fields[1]
			}
		}
	}
	runParams := genometrics.GetRunParams(testnum, mafdelta, callrate, infoscore)
	log.Printf("Params: %v\n", runParams)

	for key, value := range assaytype_filename {
		fh, err := os.Open(value)
		check(err)
		defer fh.Close()
		gr, err := gzip.NewReader(fh)
		check(err)
		defer gr.Close()
		scanner := bufio.NewScanner(gr)
		reader := bufio.NewReader(gr)
		fscanners[key] = scanner
		freaders[key] = reader
	}
	// handle file headers
	headers := make(map[string][]string)
	for assaytype, rdr := range freaders {
		headers[assaytype], _ = get_sample_headers(rdr)
		//fmt.Printf("%s hdr len = %d\n", assaytype, len(headers[assaytype]))
	}
	// Headers and combined header map
	sample_name_map, sample_posn_map := sample.MakeSamplesByAssaytype(headers)
	combocols := sample.GetCombinedSampleMap(sample_name_map)
	// combocols := sample.GetCombinedSampleMapByAssaytypes(sample_name_map, assaytype_list)
	colhdr_str, combo_names := vcfmerge.GetCombinedColumnHeaders(combocols)
	//fmt.Printf("%s\n", "combined"+"\t"+colhdr_str)
	print_headers()
	fmt.Printf("%s\n", colhdr_str)

	// read first records and capture keys (genomic positions)
	records := make(map[string][]string)
	keys := make(map[string]int64)
	varids := make(map[string]string)
	for assaytype, rdr := range freaders {
		records[assaytype], keys[assaytype], varids[assaytype] = get_next_record_slice(rdr)
	}
	var genomet genometrics.AllMetrics
	outctr := 0
	// process until all files exhausted
	for records_remain(keys) {
		records, keys, varids = check_low_key_records(records, keys, varids)
		output_from_low_key_records(records, keys, sample_posn_map, combocols, combo_names, threshold, &genomet)
		outctr += 1
		records, keys, varids = read_from_low_key_records(records, keys, freaders, varids)
	}
	log.Printf("EXIT,wrt=%d,allgeno=%d,2ol=%d,gt2ol=%d,mmc=%d,misstested=%d,missing=%d\n", outctr, genomet.AllGenoCount, genomet.TwoOverlapCount, genomet.GtTwoOverlapCount, genomet.MismatchCount, genomet.MissTestCount, genomet.MissingCount)
}

//-------------------------------------------------------------
// Get headers with column(sample) names
//-------------------------------------------------------------
func get_sample_headers(rdr *bufio.Reader) ([]string, error) {
	var sfx []string

	eof := false
	for !eof {
		text, err := rdr.ReadString('\n')
		if err == io.EOF {
			return empty_record, err
		}
		text = strings.TrimRight(text, "\n")
		if strings.HasPrefix(text, "##") {
			//fmt.Printf("%s\n", text)
		} else {
			if strings.HasPrefix(text, "#") {
				_, sfx = variant.GetVCFPrfx_Sfx(strings.Split(text, "\t"))
				break
			}
		}
	}
	return sfx, nil
}

//-------------------------------------------------------------
// Read a record from a single reader and split it to format a string slice
//-------------------------------------------------------------
func get_next_record_slice(rdr *bufio.Reader) ([]string, int64, string) {
	posn := max_posn
	data := empty_record
	varid := ""
	text, err := rdr.ReadString('\n')
	if err == nil {
		text = strings.TrimRight(text, "\n")
		data = strings.Split(text, "\t")
		posn = variant.GetPosn(data)
		varid = variant.GetVarid(data)
	}
	return data, posn, varid
}

//-------------------------------------------------------------
// If all keys are high there are no records to be read
//-------------------------------------------------------------
func records_remain(keys map[string]int64) bool {
	for _, pos := range keys {
		if pos < max_posn {
			return true
		}
	}
	return false
}

//-------------------------------------------------------------
// Cross check low key records to match alleles (when there
// are more than one)
//-------------------------------------------------------------
func check_low_key_records(records map[string][]string, keys map[string]int64,
	varids map[string]string) (map[string][]string, map[string]int64, map[string]string) {
	low_keys := get_low_keys(keys)
	key_count := 0
	if len(low_keys) > 1 {
		if key_count == 0 {
		}
	}
	return records, keys, varids
}
func output_from_low_key_records(records map[string][]string, keys map[string]int64,
	sample_posn_map map[string]map[int]string,
	combocols map[string]int, combo_names []string, threshold float64, genomet *genometrics.AllMetrics) {
	//
	low_keys := get_low_keys(keys)
	low_key_at := make([]string, 0)

	for at, _ := range low_keys {
		low_key_at = append(low_key_at, at)
	}
	sort.Strings(low_key_at)

	vcfrecords := make([][]string, 0, len(records))
	rsid := ""
	for _, at := range low_key_at {
		rsid = variant.GetVarid(records[at])
		rec := make([]string, 1, len(records[at])+1)
		rec[0] = at
		rec = append(rec, records[at]...)
		vcfrecords = append(vcfrecords, rec)
		//fmt.Printf("LOWKEY OUTPUT %s, %d\n", at, key)
	}
	var vcfd []vcfmerge.Vcfdata
	comborec := vcfmerge.Mergeslices_full(vcfrecords, vcfd, rsid, sample_posn_map, combocols, combo_names, threshold, genomet)
	rec_str := strings.Join(comborec, "\t")
	fmt.Printf("%s\n", rec_str)
}

func read_from_low_key_records(records map[string][]string, keys map[string]int64, rdrs map[string]*bufio.Reader, varids map[string]string) (map[string][]string, map[string]int64, map[string]string) {
	low_keys := get_low_keys(keys)
	for assaytype, _ := range low_keys {
		records[assaytype], keys[assaytype], varids[assaytype] = get_next_record_slice(rdrs[assaytype])
	}
	return records, keys, varids
}

func get_low_keys(keys map[string]int64) map[string]int64 {
	low_keys := make(map[string]int64)
	low_key := max_posn
	for _, pos := range keys {
		if pos < low_key {
			low_key = pos
		}
	}
	if low_key < max_posn {
		for at, pos := range keys {
			if pos == low_key {
				low_keys[at] = pos
			}
		}
	}
	return low_keys
}

func print_headers() {
	fmt.Printf("%s\n", "##fileformat=VCFv4.2")
	fmt.Printf("%s\n", "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes\">")
	fmt.Printf("%s\n", "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">")
	fmt.Printf("%s\n", "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")
	fmt.Printf("%s\n", "##INFO=<ID=RefPanelAF,Number=A,Type=Float,Description=\"Allele frequency in imputation reference panel\">")
	fmt.Printf("%s\n", "##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"Genotype dosage\">")
	fmt.Printf("%s\n", "##FORMAT=<ID=GP,Number=G,Type=Float,Description=\"Genotype posterior probabilities\">")
	fmt.Printf("%s\n", "##FORMAT=<ID=AT,Number=1,Type=String,Description=\"Assay Type\">")
	fmt.Printf("%s\n", "##INFO=<ID=TYPED,Number=0,Type=Flag,Description=\"Typed in input data\">")
}
