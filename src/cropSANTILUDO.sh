mkdir -p ./../output/

name=$1
pdfcrop ./../output/1.pdf ./../output/1.pdf 
pdfcrop ./../output/2.pdf ./../output/2.pdf
pdfcrop ./../output/3.pdf ./../output/3.pdf
pdfcrop ./../output/4.pdf ./../output/4.pdf
pdfcrop ./../output/5.pdf ./../output/5.pdf

pdfcrop ./../output/1000.pdf ./../output/1000.pdf 
pdfcrop ./../output/2000.pdf ./../output/2000.pdf
pdfcrop ./../output/3000.pdf ./../output/3000.pdf

pdfjam -q --noautoscale true --a0paper --fitpaper true --nup  3x2 ./../output/2.pdf ./../output/3.pdf ./../output/4.pdf ./../output/1000.pdf ./../output/2000.pdf ./../output/3000.pdf -o ./../output/test.pdf
pdfcrop ./../output/test.pdf ./../output/test.pdf

pdfjam -q --noautoscale true --a0paper --fitpaper true --nup  3x2 ./../output/2.pdf ./../output/3.pdf ./../output/4.pdf ./../output/1.pdf ./../output/5.pdf -o ./../output/t5.pdf
pdfcrop ./../output/t5.pdf ./../output/t5.pdf


pdfjam ./../output/t5.pdf ./../output/test.pdf --landscape -o ./../output/$name.SANTILUDO.pdf


rm -rf ./../output/1.pdf ./../output/2.pdf ./../output/3.pdf ./../output/4.pdf ./../output/5.pdf ./../output/1000.pdf ./../output/2000.pdf ./../output/3000.pdf ./../output/test.pdf ./../output/t5.pdf 
