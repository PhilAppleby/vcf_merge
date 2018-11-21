// Everything variant related
package variant

import (
	"fmt"
	"log"
	"strconv"
	"strings"
)

var sglDigitChrom map[string]int
var geno_strings = []string{"0/0", "0/1", "1/1"}

const firstGenoIdx = 9
const chrIdx = 0
const posnIdx = 1
const varIdx = 2
const refIdx = 3
const altIdx = 4
const qcIdx = 5
const filtIdx = 6
const infoIdx = 7
const fmtIdx = 8

func init() {
	sglDigitChrom = make(map[string]int)
	sglDigitChrom["atest"] = 1
	sglDigitChrom["atest2"] = 1
	sglDigitChrom["atest3"] = 1
	sglDigitChrom["atest4"] = 1
	sglDigitChrom["atest5"] = 1
	sglDigitChrom["affy"] = 1
	sglDigitChrom["illumina"] = 1
	sglDigitChrom["broad"] = 1
}

func check(e error) {
	//  fmt.Println("check")
	if e != nil {
		fmt.Println("err != nil")
		log.Fatal(e)
	}
}

func AppendToFmt(prfx []string, add_str string) []string {
	prfx[fmtIdx] = prfx[fmtIdx] + ":" + add_str
	return prfx
}

func NormaliseChromosome(prfx []string) []string {
	if prfx[chrIdx][0] == '0' {
		prfx[chrIdx] = prfx[chrIdx][1:]
	}
	return prfx
}

//------------------------------------------------------------------------------
// get geno based on threshold
//------------------------------------------------------------------------------
func Get_geno(geno string, threshold float64, probidx int) string {

	mprob, max_prob_idx, genoarray := MaxProb(geno, probidx)
	//fmt.Printf("GGENO %s,%f,%f,%d\n", geno, mprob, threshold, probidx)
	if mprob < threshold {
		genoarray[0] = "./."
	} else {
		genoarray[0] = geno_strings[max_prob_idx]
	}
	return strings.Join(genoarray, ":")
}

//------------------------------------------------------------------------------
// maxprob test for a genotype
//------------------------------------------------------------------------------
func MaxProb(geno string, probidx int) (float64, int, []string) {
	g := strings.Split(geno, ":")
	//if len(g) < (probidx + 1) {
	//	fmt.Printf("GENO IDX PROBLEM %d, %s\n", probidx, geno)
	//	return 0.0, -9, g
	//}
	probs := strings.Split(g[probidx], ",")
	max_prob := 0.0
	max_prob_idx := -9

	for i, prob := range probs {
		probf, _ := strconv.ParseFloat(prob, 64)
		if probf > max_prob {
			max_prob = probf
			max_prob_idx = i
		}
	}
	return max_prob, max_prob_idx, g
}

//------------------------------------------------------------------------------
// ---- a group of convenience fns to get array fields
//------------------------------------------------------------------------------
func GetVCFPrfx_Sfx(recslice []string) ([]string, []string) {
	return recslice[:firstGenoIdx], recslice[firstGenoIdx:]
}

func GetVarid(recslice []string) string {
	return recslice[varIdx]
}

func GetPosn(recslice []string) int64 {
	value, _ := strconv.ParseInt(recslice[posnIdx], 0, 64)
	return value
}
func GetAlleles(recslice []string) (string, string) {
	return recslice[refIdx], recslice[altIdx]
}

func GetProbidx(recslice []string) int {
	return getStrIdx(recslice[fmtIdx], "GP")
}

func getStrIdx(str string, match_str string) int {
	str_arr := strings.Split(str, ":")
	for i, mstr := range str_arr {
		if mstr == match_str {
			return i
		}
	}
	return -9
}
