# IBDMUT

**Reference:** P.F. Palamara, L. Francioli, P. Wilton, G. Genovese, A. Gusev, H. Finucane, S. Sankararaman, GoNL Consortium, S. Sunyaev, P. de Bakker, J. Wakeley, I. Pe’er, A. Price, "[Leveraging Distant Relatedness to Quantify Human Mutation and Gene-Conversion Rates](www.cell.com/ajhg/abstract/S0002-9297(15)00408-5)", American Journal of Human Genetics, 2015.

**Contact:** Pier Palamara, ppalama AT hsph DOT harvard DOTAGAIN edu

IBDMUT Infers mutation and gene conversion rates using identical-by-descent segments

### Usage

Options:

    --demography [file] (demographic model with format "generation	size" for each line. Sizes are haploid)
    --plink [fileRoot] (plink file root; expects .ped.gz, .map, and .frq)
    --match [file] (match.gz file in GERMLINE format)
	
Some of the optional or alternative flags:

    --MaAFRegression [intValueFromCount] [intValueStep] [intValueToCount] (from MAF counts, interval MAF counts, to MAF counts)
    --offsetCM [value] (distance to be excluded from edges; default: 0.0)
    --lenRange [fromLen] [toLen] (default: 2.0 5.0)
    --plinkList [file] (substitutes --plink and --match for many files)
    --jackknife (only if --plinkList is used with several independent regions)
    --saveBin [suffix] (saves a binary file, which will load much faster than the ped.gz file)
    --loadBin [suffix] (load a binary file)
    --mask [file] (bed file with regions to be included in analysis)
    --maskList [file] (substitutes --mask for many files)
    --threads [value] (default: 1)

### Example

    java -jar IBDMUT.jar \
    	--plinkList EXAMPLES/plinkList.txt \
    	--demography EXAMPLES/10K.demo \
    	--lenRange 1.6 5.0 \
    	--offsetCM 0.5 \
    	--MaAFRegression 50 5 250 \
    	--saveBin .save \
    	--jackknife \
    	--threads 4

### Dependencies

(all in lib folder)

Kryo: https://github.com/EsotericSoftware/kryo

Minlog: https://github.com/EsotericSoftware/minlog/

Objenesis: https://github.com/easymock/objenesis

Reflactasm: https://github.com/EsotericSoftware/reflectasm
