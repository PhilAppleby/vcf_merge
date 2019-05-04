package genometrics

import (
	//	"fmt"
	"strconv"
	"strings"
	"variant"
)

type AllMetrics struct {
	AllGenoCount      int
	UniqueGenoCount   int
	OverlapTestCount  int
	TwoOverlapCount   int
	GtTwoOverlapCount int
	MismatchCount     int
	MissTestCount     int
	MissingCount      int
}

type RunParameters struct {
	TestNum   int
	MafDelta  float64
	CallRate  float64
	InfoScore float64
}

// caller passes a string array representing a whole VCF
// record, including prefix
func Hwe_exact_for_record(rec []string, threshold float64) float64 {
	homref, homalt, het, _, _, _, _ := get_genotype_counts(rec, threshold)
	return SNPHWE(het, homref, homalt)
}

// return all SNP metrics
// CR, RAF, AAF, MAF, HWE_P
func Metrics_for_record(rec []string, threshold float64) (float64, float64,
	float64, float64, float64, int, int, int, int, int, int, float64) {
	homref, homalt, het, n, miss, dot, refPAF := get_genotype_counts(rec, threshold)
	cr := float64(homref+het+homalt) / float64(n)
	raf := float64(2*homref+het) / float64(2*n)
	aaf := float64(2*homalt+het) / float64(2*n)
	maf := aaf
	if raf < aaf {
		maf = raf
	}
	obs_homc := homref
	obs_homr := homalt
	if homalt > homref {
		obs_homc = homalt
		obs_homr = homref
	}
	return cr, raf, aaf, maf, SNPHWE(het, homref, homalt), het, obs_homc, obs_homr, n, miss, dot, refPAF
}

func GetRunParams(testnum string, mafdelta string, callrate string, infoscore string) RunParameters {
	var runParams RunParameters

	runParams.TestNum = 0
	runParams.MafDelta = 0.0
	runParams.CallRate = 0.0
	runParams.InfoScore = 0.0

	if testnum != "" {
		runParams.TestNum, _ = strconv.Atoi(testnum)
	}
	if mafdelta != "" {
		runParams.MafDelta, _ = strconv.ParseFloat(mafdelta, 32)
	}
	if callrate != "" {
		runParams.CallRate, _ = strconv.ParseFloat(callrate, 32)
	}
	if infoscore != "" {
		runParams.InfoScore, _ = strconv.ParseFloat(infoscore, 32)
	}

	return runParams
}

// Following is translated from the 'C' routine from UMICH
// Original comments:
//
// "This code implements an exact SNP test of Hardy-Weinberg Equilibrium as described in
// Wigginton, JE, Cutler, DJ, and Abecasis, GR (2005) A Note on Exact Tests of
// Hardy-Weinberg Equilibrium. American Journal of Human Genetics. 76: 000 - 000
//
// Written by Jan Wigginton"
//
func SNPHWE(obs_hets int, obs_hom1 int, obs_hom2 int) float64 {
	obs_homc := obs_hom1
	obs_homr := obs_hom2

	if obs_hom2 > obs_hom1 {
		obs_homc = obs_hom2
		obs_homr = obs_hom1
	}

	rare_copies := 2*obs_homr + obs_hets
	genotypes := obs_hets + obs_homc + obs_homr

	het_probs := make([]float64, rare_copies+1, rare_copies+1)

	// start at the midpoint
	mid := rare_copies * (2*genotypes - rare_copies) / (2 * genotypes)
	//fmt.Printf("mid %d, rare %d\n", mid, rare_copies)

	if ((rare_copies & 1) ^ (mid & 1)) != 0 {
		mid += 1
		//fmt.Printf("add 1 to mid %d\n", mid)
	}

	curr_homr := (rare_copies - mid) / 2
	curr_homc := genotypes - mid - curr_homr

	het_probs[mid] = 1.0
	sum := 1.0

	for curr_hets := mid; curr_hets > 1; curr_hets -= 2 {
		het_probs[curr_hets-2] = het_probs[curr_hets] * float64(curr_hets) * float64(curr_hets-1.0) / (4.0 * float64(curr_homr+1.0) * float64(curr_homc+1.0))
		sum += het_probs[curr_hets-2]
		// 2 fewer heterozygotes for next iteration -> add one rare, one common homozygote
		curr_homr += 1
		curr_homc += 1
	}

	curr_homr = (rare_copies - mid) / 2
	curr_homc = genotypes - mid - curr_homr

	for curr_hets := mid; curr_hets <= rare_copies-2; curr_hets += 2 {
		het_probs[curr_hets+2] = het_probs[curr_hets] * 4.0 * float64(curr_homr) * float64(curr_homc) / (float64(curr_hets+2.0) * float64(curr_hets+1.0))
		sum += het_probs[curr_hets+2]
		// add 2 heterozygotes for next iteration -> subtract one rare, one common homozygote
		curr_homr -= 1
		curr_homc -= 1
	}

	for i, _ := range het_probs {
		het_probs[i] /= sum
	}

	p_hwe := 0.0

	for i := 0; i <= rare_copies; i++ {
		if het_probs[i] <= het_probs[obs_hets] {
			p_hwe += het_probs[i]
		}
	}
	if p_hwe > 1.0 {
		p_hwe = 1.0
	}
	return p_hwe

}

func get_genotype_counts(rec []string, threshold float64) (int, int, int, int, int, int, float64) {
	homr := 0
	homa := 0
	het := 0
	n := 0
	miss := 0
	dot := 0

	prfx, sfx := variant.GetVCFPrfx_Sfx(rec)
	probidx := variant.GetProbidx(prfx)
	refPAF := variant.GetRefPanelAf(prfx)

	for _, geno := range sfx {
		if geno != "." {
			n += 1
			geno = variant.Get_geno(geno, threshold, probidx)
			geno_a := strings.Split(geno, ":")
			if geno_a[0] == "0/0" {
				homr += 1
			}
			if geno_a[0] == "0/1" {
				het += 1
			}
			if geno_a[0] == "1/0" {
				het += 1
			}
			if geno_a[0] == "1/1" {
				homa += 1
			}
			if geno_a[0] == "./." {
				miss += 1
			}
		} else {
			dot += 1
		}
	}

	return homr, homa, het, n, miss, dot, refPAF
}
