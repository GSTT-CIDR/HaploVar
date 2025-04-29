#!/bin/bash


echo "do you want phylogenies? (yes/no; default = no)"
read -t 5 treeresponse
if [ -z "$treeresponse" ]; then
    treeresponse="no"  
fi

function setlogs {
    # Give sudo permissions
    sudo ls

    # Define the log file
    log_file="haplotype_script_log.txt"

    # Redirect stdout to the terminal and the log file
    exec > >(tee -a "$log_file")

    # Redirect stderr to stdout (so that errors are also logged)
    exec 2>&1
}

function extractrunname {
    # Find the run name by extracting filename prior to BC
    for file in *.fastq; do
        runname="${file%%BC*}"
    done
}

function extractq30intersect {
    # Process FASTQ files
    for fastq in *fastq; do
        cat "$fastq" | NanoFilt -q 30 | minimap2 -t 16 -ax map-ont SARS-CoV-2.reference.fasta - -t 16 | samtools sort - | bedtools intersect -wo -a spike.bed -b - | awk '$5 < 21564' - | awk '$6 > 25384' - | cut -f7 - | seqtk subseq "$fastq" - > "${fastq//.fastq}.q30.intersect.fastq"
    done
}

function createbam {
    # Map q30 reads to BAM
    for i in *.q30.intersect.fastq; do
        minimap2 -t 16 -ax map-ont SARS-CoV-2.reference.fasta "$i" -t 16 | samtools sort - > "${i//.fastq}.bam"
    done
    for i in *bam; do
        samtools index "$i"
    done
}

function logq30spanningreads {
    # Check how many q30 intersecting reads and log them
    for i in *.q30.intersect.fastq; do
        echo -n "$i " >> q30_spanning_reads.txt && cat "$i" | wc -l | awk '{print $1/4}' >> q30_spanning_reads.txt
    done
}

function runfreebayes {
    # Run FreeBayes to find variants
    
    for i in *bam; do
        freebayes -f SARS-CoV-2.reference.fasta -F 0.01 -C 20 --pooled-continuous "$i" > "${i//.bam}freebayes_0.01.tsv"
    done
    
}

function findpositionsofvariance {    
 #From the FreeBayes output TSV, pull out non-fixed positions of variance with qual > 100
    for i in *tsv; do
        awk -F'\t' '$10 !~ /1\/1/ && $6 >= 100 {chars_in_d = length($4) - 1; for (i = $2; i <= $2 + chars_in_d; i++) print i}' "$i" > "${i//.tsv}_positionofvariance.txt"
    done
}

function housekeeping {
    # Tidy the output into sample specific folders
    for i in {01..24}; do
        if compgen -G "*BC$i*" > /dev/null; then
            mkdir -p "$runname"BC"$i"
            mv *BC$i* "$runname"BC"$i"/
        fi
    done
}

function haplotyping {
    # Enter each directory, create fasta from fastq for haplotype processing
    for dir in "${runname}"BC{01..24}; do
        if [ -d "$dir" ]; then
            echo "Entering directory $dir"
            cd "$dir" || exit
            echo -e "looking for positions of variance file (freebayes 0.01)"
            positions_file=$(ls *_positionofvariance.txt | head -n 1)
            echo "Got it"
            outputfolder=haplotype_output
            echo "check fastq"
            fastq=(*duplex.q30.intersect.fastq)
            echo -e "set output file"
            output_file="${fastq%.duplex.q30.intersect.fastq}"
            echo "Got it"
            echo -e "\nGo"
            # Check if the positions file is empty
            if [[ ! -s "$positions_file" ]] && [[ $(cat "$fastq" | wc -l | awk '{print $1/4}') -ge 10 ]]; then
                echo "${fastq//.duplex.q30.intersect.fastq}: the positions file is empty; single haplotype. Run NextFlow" > "${fastq//.duplex.q30.intersect.fastq}.log.txt"
                conda activate nextflow
                sudo nextflow run wf-artic/main.nf --scheme_name SARS-CoV-2 --scheme_version Midnight-IDT/V1 --fastq "$fastq" --sample "${fastq//.fastq}" --max_len 4300 --out_dir "${fastq//.fastq}_nextflow" --basecaller_cfg dna_r10.4.1_e8.2_400bps_sup@v3.5.2
                conda deactivate
                cat *nextflow/*consensus.fasta > "${fastq//.duplex.q30.intersect.fastq}_haplotype_consensus.fasta"
                cd ../
            elif [[ ! -s "$positions_file" ]] && [[ $(cat "$fastq" | wc -l | awk '{print $1/4}') -lt 10 ]]; then
            	echo "${fastq//.duplex.q30.intersect.fastq}: fails; not enough (<10) Q30 spanning reads" > "${fastq//.duplex.q30.intersect.fastq}.log.txt"
            	cd ../
            else
            	echo "${fastq//duplex.q30.intersect.fastq}: multiple haplotypes found. Start haplotyping" > "${fastq//.duplex.q30.intersect.fastq}.log.txt"
                source .bashrc
                cat "$fastq" | NanoFilt -q 30 | minimap2 -t 16 -ax map-ont SARS-CoV-2.reference.fasta - -t 16 | samtools sort - | bedtools intersect -wo -a spike.bed -b - | awk '$5 < 21564' - | awk '$6 > 25384' - | cut -f7 - | seqtk subseq "$fastq" - | seqtk seq -a - | minimap2 -t 16 -ax map-ont SARS-CoV-2.reference.fasta - -t 16 | gofasta sam toMultiAlign --reference SARS-CoV-2.reference.fasta - > "${fastq//.fastq/}.intersect.fasta"
                cat "${fastq//.fastq/}.intersect.fasta" | minimap2 -t 16 -a SARS-CoV-2.reference.fasta -  -t 16 | gofasta sam toMultiAlign --reference SARS-CoV-2.reference.fasta - | awk '{print substr($1,1,74)}{print substr($0,21563,3822)}'  |  awk '!/^$/' | sed '2~3d' > "${fastq//.fastq/}.spikeonly.fasta"
                echo "Extract positions of variance"
                {
                    while read -r line; do
                        if [[ $line =~ ^\> ]]; then
                            echo "$line"
                        else
                            sequence="$line"
                            while read -r pos; do
                                character=$(echo "$sequence" | awk -v pos="$pos" '{print substr($0, pos, 1)}')
                                echo -n "$character"
                            done < "$positions_file"
                            echo
                        fi
                    done < "${fastq//.fastq/}.intersect.fasta"
                } > "$output_file"
                echo "Positions of variance captured"
                echo "Convert fasta to tsv for further processing"
                positions_of_variance_fasta="$output_file"
                output_tsv="$positions_of_variance_fasta"_positionsofinterest.tsv
                header=""
                sequence=""
                while IFS= read -r line; do
                    if [[ $line =~ ^\> ]]; then
                        if [[ -n $header ]]; then
                            echo -e "$header\t$sequence" >> "$output_tsv"
                        fi
                        header="$line"
                        sequence=""
                    else
                        sequence="$sequence$line"
                    fi
                done < "$positions_of_variance_fasta"
                if [[ -n $header ]]; then
                    echo -e "$header\t$sequence" >> "$output_tsv"
                fi
                echo "Conversion completed. Output saved to $output_tsv"
                cut -f 2 "$output_tsv" | sort | uniq > "${output_tsv//.tsv}"_uniques.tsv
                echo "Unique haplotypes identified"
                echo "Counting occurrences of each haplotype"
                grep -wf <(cut -f1 *uniques.tsv) "$output_tsv" | cut -f2 | sort | uniq -c | sort -nrk1,1 | awk -v OFS='\t' '{ print $2, $1 }' > "${output_tsv//.tsv}"_counts.tsv
                echo "split reads into haplotype groups"
                # Counter for row number
		rowNumber=0
		# Read each line of the TSV file
		while IFS=$'\t' read -r searchString frequency; do
    			# Increment row number
    			((rowNumber++))
    
    			echo -e "${searchString}\t${frequency}\t${rowNumber}" >> haplotyperank.txt

    			# Search searchString in the second column of datafile
    			grep -F -- "$searchString" *positionsofinterest.tsv | while IFS=$'\t' read -r dataValue dataKey; do
        			# Check if the searchString matches the dataKey
        			if [ "$searchString" = "$dataKey" ]; then
            				# Append the dataValue to a file named after the row number
            				echo "$dataValue" >> "haplotype_${rowNumber}.txt"
        			fi
    			done
		done < *counts.tsv
                perl -i -pe 's/>//g' haplotype_*.txt
                echo "done"
                echo "subseq the fastq into haplotype groups"
                for haplotypetext in haplotype_*.txt; do 
                	rank="${haplotypetext#haplotype_}"
    			rank="${rank%.txt}"
    			linecount=$(wc -l < "$haplotypetext" | xargs)
                	seqtk subseq "$fastq" "$haplotypetext" > "${fastq//.fastq}.haplotype.$rank.$linecount.fastq"; 
                done
                echo "subseq done"
                echo "tidy output"
                mkdir -p $outputfolder
                mkdir -p haplotype_working
                mv *tsv haplotype_working/
                mv $outputfolder/*tsv haplotype_working/
                mv *.haplotype.*.fastq $outputfolder
                cd $outputfolder
                echo "done"
                echo "run nextflow and artic standalone"
                sleep 5
                for i in *.fastq; do
                    source conda.sh
                    read_count=$(( $(cat "$i" | wc -l) / 4 | bc))
                    if (($read_count > 19)); then
                        echo -e "running nextflow for haplotypes >= 20 reads"
                        sudo nextflow run wf-artic/main.nf --scheme_name SARS-CoV-2 --scheme_version Midnight-IDT/V1 --fastq "$i" --sample "${i//.fastq}" --artic_threads 16 --max_len 4300 --out_dir "${i//.fastq}_nextflow" --basecaller_cfg dna_r10.4.1_e8.2_400bps_sup@v3.5.2                     
                    elif (($read_count > 9 && $read_count < 20)); then
                        read_count=$(( $(cat "$i" | wc -l) / 4 | bc))
                        echo -e "running artic standalone for haplotypes >=10 & <20"
                        conda activate artic
                        mkdir -p artic_standalone
                        cd artic_standalone
                        artic minion --medaka --normalise 200 --threads 16 --scheme-directory SARS-CoV-2 --read-file  ../"$i"    Midnight-IDT/V1 "${i//.fastq}"  --medaka-model r1041_e82_400bps_sup_variant_g615
                        conda deactivate
                        cd ../
                    fi
                done
                echo "tidying"
                mkdir -p haplotype_fastqs
                mv *fastq  haplotype_fastqs
                echo "combine haplotype consensuses"
                source .bashrc
                cat *nextflow/all_consensus.fasta artic_standalone/*.consensus.fasta > "${fastq//.duplex.q30.intersect.fastq}_haplotype_consensus.fasta"
                minimap2 -t 16 -a SARS-CoV-2.reference.fasta "${fastq//.duplex.q30.intersect.fastq}_haplotype_consensus.fasta" -t 16 | gofasta sam toMultiAlign --reference SARS-CoV-2.reference.fasta - | awk '{print substr($1,1,74)}{print substr($0,21563,3822)}' -  | awk '!/^$/' | sed '2~3d' > "${fastq//.duplex.q30.intersect.fastq}_haplotype_consensus.spikeonly.fasta"
                echo "exit"
                cd ../../
            fi
    	fi
    done
}
function phylogeny {
	for dir in ${runname}BC{01..24}; do
        	if [ -d "$dir" ]; then
        		echo "Entering directory $dir"
        		cd "$dir" || exit
        		intersectforphylogeny=(*duplex.q30.intersect.spikeonly.fasta)
			iqtree2 -s $intersectforphylogeny -m K3Pu+F+I+G4  -alrt 1000 -T AUTO
			cd ../
		fi
	done	 
}


function annotatephylogeny {
	source .bashrc
	for dir in "${runname}"BC{01..24}; do
        	if [ -d "$dir" ]; then
			samplename=$(basename $dir)
			cd "$dir/haplotype_output/haplotype_fastqs/"
			# Write header to the output tsv
			echo "start phylogeny annotation for $samplename"
			echo -e "taxa\tcluster\tsample\tsample_cluster" > "$samplename.annotation.tsv"
			for i in *fastq; do
				echo "find fastqs have more than 10 reads"
				if [[ $(cat "$i" | wc -l | awk '{print $1/4}') -ge 10 ]]; then	
					# Iterate over all fastq files in the current directory
					echo "start finding the read ids"
					# Extract read IDs and count them
    					awk 'NR%4==1' "$i" | while read id; do
        				# Count the number of read IDs in the fastq file
        				count=$(awk 'NR%4==1' "$i" | wc -l)
        				echo -e "$id\t$count" >> "$samplename.annotation.tsv_temp.tsv"
    			
					sed -i 's/@//g' "$samplename.annotation.tsv_temp.tsv"
					sed -i 's/;/_/g' "$samplename.annotation.tsv_temp.tsv"

					awk -F'\t' -v OFS='\t' -v token="$samplename" '{print $0, token}' "$samplename.annotation.tsv_temp.tsv" >> "$samplename.annotation.tsv_temp2.tsv"

					awk -F'\t' 'BEGIN {OFS="\t"} {print $0, $3 "_" $2}' "$samplename.annotation.tsv_temp2.tsv" >> "$samplename.annotation.tsv"


					rm "$samplename.annotation.tsv_temp.tsv"
					rm "$samplename.annotation.tsv_temp2.tsv"
					done
				fi
			done
			mv "$samplename.annotation.tsv" ../../
			cd ../../../
		fi
		
	done	
}



# Call the functions sequentially
setlogs
extractrunname
extractq30intersect
createbam
logq30spanningreads
runfreebayes
findpositionsofvariance
housekeeping
haplotyping
if [ $treeresponse = "yes" ]; then phylogeny; fi
if [ $treeresponse = "yes" ]; then annotatephylogeny; fi


