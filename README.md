# IBDMUT
Infer mutation and gene conversion rates using identical-by-descent segments

Contact: Pier Palamara, ppalama AT hsph DOT harvard DOTAGAIN edu

### Usage

Options:

    --demography [file] (demographic model with format "generation	size" for each line)
    --plink [fileRoot] (plink file root; expects .ped.gz, .map, and .frq)
    --match [file] (match.gz file in GERMLINE format)
	
Some of the optional or alternative flags:

    --plinkList [file] (substitutes --plink and --match for many files)
    --maskList [file] (substitutes --mask for many files)
    --lenRange [fromLen] [toLen] (default: 2.0 5.0)
    --jackknife (only if --plinkList is used with several independent regions)
    --saveBin [suffix] (saves a binary file, which will load much faster than the ped.gz file)
    --loadBin [suffix] (load a binary file)
    --offsetCM [value] (distance to be excluded from edges; default: 0.0)
    --mask [file] (bed file with regions to be included in analysis)
    --threads [value] (default: 1)
    --MaAFRegression [intValueFromCount] [intValueStep] [intValueToCount] (from MAF counts, interval MAF counts, to MAF counts)

### Example

    java -jar IBDMUT.jar \
    	--plinkList EXAMPLES/plinkList.txt \
    	--demography EXAMPLES/10K.demo \
    	--lenRange 1.6 5.0 --offsetCM 0.5 \
    	--MaAFRegression 50 5 250 \
    	--jackknife \
    	--threads 4

### Dependencies

(all in lib folder)

Kryo: https://github.com/EsotericSoftware/kryo

Minlog: https://github.com/EsotericSoftware/minlog/

Objenesis: https://github.com/easymock/objenesis

Reflactasm: https://github.com/EsotericSoftware/reflectasm

### Source code

included in jar file

### Citation

- P.F. Palamara et al. "Leveraging distant relatedness to quantify human mutation and gene conversion rates", American Journal of Human Genetics, 2015.
