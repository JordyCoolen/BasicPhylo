
FILES="/ifs/home/jordy/paul/pipelines/main/data/new_fasta/*"
for f in $FILES
do
	fullname="${f##*/}" 
	filename="${fullname%.*}"
	ska fasta -k 15 -o data/temp/$filename data/new_fasta/$fullname

done
