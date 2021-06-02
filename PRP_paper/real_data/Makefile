all: 	 
	if [ ! -d "output" ]; then mkdir output; fi
	Rscript scripts/analyze_RPP.R
	Rscript scripts/analyze_CAD.R > output/CAD_rst.txt
clean:
	rm -rf output/

