
 
docker run -v $(pwd):/research cryanking/temperature_container  R --slave -e 'options(width=2000); writeLines(toLatex(sessionInfo()),"rinfo.tex") ; write.csv(installed.packages()[, c("Package", "Version", "Built", "LibPath")], "/research/r_packages.txt")'

echo $(date +'%m/%d/%Y %r') > build.txt
echo -n "build commit: " >> build.txt
git rev-parse --verify HEAD >> build.txt

echo "input md5:" >> build.txt
md5sum -- ob_temp_working_data.csv >> build.txt 

