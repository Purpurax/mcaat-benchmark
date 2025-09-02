#!/bin/bash

usage() {
    echo "Usage: $0 <mcaat_binary> <genome_crispr_combination.csv> [OPTIONS]"
    echo ""
    echo "Requirements:"
    echo "  unzip     Installed through apt-get install unzip"
    echo "  datasets  See more on https://github.com/ncbi/datasets"
    echo "  iss       See more on https://github.com/HadrienG/InSilicoSeq"
    echo ""
    echo "Required arguments:"
    echo "  mcaat_binary                   The compiled application https://github.com/RNABioInfo/mcaat"
    echo "  genome_crispr_combination.csv  Extracted from sql dump, see instructions on https://github.com/RNABioInfo/mcaat/todo"
    echo "  --ids <id1> or more ids"
    echo ""
    echo "Optional arguments:"
    echo "  --ram <size>            Memory limit (e.g., 4G, 500M). Default: 80% of system RAM"
    echo "  --threads <count>       Number of threads. Default: 2"
    echo "  --output_folder <path>  Output directory for results"
    echo "  --quick                 Benchmark using a reduced genome"
    echo "  --ids <n_samples> || <id1> [id2 ...]   Space-separated IDs for benchmarking, or a number of random samples"
    echo ""
}

verify_requirements() {
    if ! command -v unzip >/dev/null 2>&1; then
        echo "Error: 'unzip' not found, please install over apt-get install unzip"
        exit 1
    fi

    if ! command -v datasets >/dev/null 2>&1; then
        echo "Error: 'datasets' not found, please install on https://github.com/ncbi/datasets"
        exit 1
    fi

    if ! command -v iss >/dev/null 2>&1; then
        echo "Error: 'iss' not found, please install on https://github.com/HadrienG/InSilicoSeq"
        exit 1
    fi
}

verify_arguments() {
    if [ "$#" -lt 4 ]; then
        usage
        exit 1
    fi

    MCAAT_BINARY="$1"
    GENOME_CRISPR_CSV="$2"

    if [ ! -f "$MCAAT_BINARY" ]; then
        echo "Error: mcaat binary file '$MCAAT_BINARY' not found"
        exit 1
    else
        echo "  ▸ File mcaat binary was found"
    fi

    if [ ! -f "$GENOME_CRISPR_CSV" ]; then
        echo "Error: genome_crispr_combination.csv file '$GENOME_CRISPR_CSV' not found"
        exit 1
    else
        echo "  ▸ File genome_crispr_combination.csv was found"
    fi

    RAM=$(free -m | awk 'NR==2{printf "%.0fM", $7*0.8}')
    THREADS=2
    OUTPUT_FOLDER="benchmark"
    QUICK=false
    IDS=""

    shift 2
    if [ "$#" -gt 0 ]; then
        echo "  ▸ Parsing optional arguments:"
    fi
    while [ "$#" -gt 0 ]; do
        case "$1" in
            --ram)
                if [ -z "$2" ]; then
                    echo "Error: --ram requires a value"
                    exit 1
                fi
                RAM="$2"
                echo "    ▸ RAM is set to $RAM"
                shift 2
                ;;
            --threads)
                if [ -z "$2" ]; then
                    echo "Error: --threads requires a value"
                    exit 1
                fi
                THREADS="$2"
                echo "    ▸ Threads are set to $THREADS"
                shift 2
                ;;
            --output_folder)
                if [ -z "$2" ]; then
                    echo "Error: --output_folder requires a value"
                    exit 1
                fi
                OUTPUT_FOLDER="$2"
                echo "    ▸ Outputting data in $OUTPUT_FOLDER"
                shift 2
                ;;
            --quick)
                QUICK=true
                echo "    ▸ Quick mode enabled"
                shift
                ;;
            --ids)
                shift
                while [ "$#" -gt 0 ] && [[ "$1" != --* ]]; do
                    IDS="$IDS $1"
                    shift
                done
                echo "    ▸ Benchmarking ids $IDS"
                ;;
            *)
                echo "Error: Unknown argument '$1'"
                usage
                exit 1
                ;;
        esac
    done

    if [ -z "$IDS" ]; then
        echo "Error: No ID specified. Please provide at least one ID using --ids."
        exit 1
    fi
}

setting_up_environment() {
    GENOMES_FOLDER="${OUTPUT_FOLDER}/genomes"
    READS_FOLDER="${OUTPUT_FOLDER}/reads"
    EXPECTED_CRISPRS_FOLDER="${OUTPUT_FOLDER}/crispr_sequences"
    RESULTS_FOLDER="${OUTPUT_FOLDER}/results"

    for folder in "$OUTPUT_FOLDER" "$GENOMES_FOLDER" "$READS_FOLDER" "$EXPECTED_CRISPRS_FOLDER" "$RESULTS_FOLDER"; do
        if [ ! -d "$folder" ]; then
            mkdir -p "$folder"
            if [ $? -eq 0 ]; then
                echo "  ▸ Created folder: $folder"
            else
                echo "Error: Failed to create folder '$folder'"
                exit 1
            fi
        else
            echo "  ▸ Folder exists: $folder"
        fi
    done

    if echo "$IDS" | grep -Eq '^[[:space:]]*[0-9]+[[:space:]]*$'; then
        NUM_SAMPLES=$(echo "$IDS" | tr -d '[:space:]')
        echo "  ▸ Converting numeric IDS to random samples"
        
        local total_lines=$(wc -l < "$GENOME_CRISPR_CSV")
        local available_samples=$((total_lines - 2))
        
        if [ "$NUM_SAMPLES" -ge "$available_samples" ]; then
            IDS=$(tail -n +2 "$GENOME_CRISPR_CSV" | awk -F',' 'NF > 0 && $1 != "" {print $1}' | tr '\n' ' ')
            NUM_SAMPLES=$available_samples
        else
            IDS=""
            local selected_ids=()
            local attempts=0
            local max_attempts=$((NUM_SAMPLES * 10))
            
            while [ ${#selected_ids[@]} -lt $NUM_SAMPLES ] && [ $attempts -lt $max_attempts ]; do
                local random_line_num=$((RANDOM % (total_lines - 1) + 2))
                local sample_id=$(sed -n "${random_line_num}p" "$GENOME_CRISPR_CSV" | awk -F',' '{print $1}' | xargs)
                
                local found=false
                for selected in "${selected_ids[@]}"; do
                    if [ "$selected" = "$sample_id" ]; then
                        found=true
                        break
                    fi
                done
                
                if [ "$found" = false ]; then
                    selected_ids+=("$sample_id")
                else
                    local current_line=$random_line_num
                    while true; do
                        current_line=$((current_line + 1))
                        if [ $current_line -gt $total_lines ]; then
                            current_line=2
                        fi
                        
                        sample_id=$(sed -n "${current_line}p" "$GENOME_CRISPR_CSV" | awk -F',' '{print $1}' | xargs)
                        
                        found=false
                        for selected in "${selected_ids[@]}"; do
                            if [ "$selected" = "$sample_id" ]; then
                                found=true
                                break
                            fi
                        done
                        
                        if [ "$found" = false ]; then
                            selected_ids+=("$sample_id")
                            break
                        fi
                        
                        if [ $current_line -eq $random_line_num ]; then
                            break
                        fi
                    done
                fi
                attempts=$((attempts + 1))
            done
            
            IDS=$(printf "%s " "${selected_ids[@]}")
        fi
        
        IDS=$(echo "$IDS" | xargs)
        echo "  ▸ Selected random IDs: $IDS"
    else
        NUM_SAMPLES=$(( $(echo "$IDS" | tr -cd ' ' | wc -c) + 1 ))
    fi
    echo "  ▸ Running through $NUM_SAMPLES samples"
}

get_random_id_sample() {
    if [ ! -f "$GENOME_CRISPR_CSV" ]; then
        echo "Error: Cannot open CSV file '$GENOME_CRISPR_CSV'"
        return 1
    fi
    
    local total_lines=$(wc -l < "$GENOME_CRISPR_CSV")
    if [ "$total_lines" -le 2 ]; then
        echo "Error: CSV file must have more than 2 lines"
        return 1
    fi
    
    local middle_lines=$((total_lines - 2))
    local random_line_num=$((RANDOM % middle_lines + 2))
    local sample_line=$(sed -n "${random_line_num}p" "$GENOME_CRISPR_CSV")
    
    echo "$sample_line" | awk -F',' '{print $1}' | xargs
}

download_genomes() {
    local success_count=0
    local total_count=0
    local successful_ids=""
    local pids=()
    local id_array=($BUCKET_IDS)

    for id in $BUCKET_IDS; do
        if [ -n "$id" ]; then
            ((total_count++))
            (
                if download_single_genome "$id" 2>/dev/null; then
                    exit 0
                else
                    exit 1
                fi
            ) &
            pids+=($!)
        fi
    done

    local i=0
    for pid in "${pids[@]}"; do
        wait "$pid"
        if [ $? -eq 0 ]; then
            ((success_count++))
            successful_ids="$successful_ids ${id_array[$i]}"
        fi
        ((i++))
    done

    BUCKET_IDS=$(echo "$successful_ids" | xargs)

    if [ "$success_count" -ne 0 ]; then
        echo "    ▸ ✔️  Downloaded ${success_count}/${total_count} genomes"
        return 0
    else
        echo "    ▸ ❌  Downloaded ${success_count}/${total_count} genomes"
        return 1
    fi
}

download_single_genome() {
    local id="$1"
    local genome_file="${GENOMES_FOLDER}/${id}.fasta"
    local genome_zip="${GENOMES_FOLDER}/${id}.zip"
    local genome_tmp_unzipped_folder="${GENOMES_FOLDER}/${id}_unzipped"

    local genbank
    genbank=$(extract_genbank_from_csv "$id")
    if [ $? -ne 0 ]; then
        echo "Unable to find the genbank id of the genome"
        return 1
    fi

    local genome_description
    genome_description=$(extract_description_from_csv "$id")
    if [ $? -ne 0 ]; then
        echo "Unable to find the description of the genome"
        return 1
    fi

    # Download full fasta file
    local fasta_file
    if datasets download genome accession "$genbank" --include genome --filename $genome_zip >/dev/null 2>&1; then
        if unzip -q "${genome_zip}" -d "$genome_tmp_unzipped_folder" >/dev/null 2>&1; then
            rm -f "${genome_zip}"

            fasta_file=$(find "$genome_tmp_unzipped_folder" -name "*.fna" -o -name "*.fasta" -o -name "*.fa" | head -n 1)

            if [ -z "$fasta_file" ] || [ ! -f "$fasta_file" ]; then
                rm -rf "$genome_tmp_unzipped_folder"

                echo "The zip doesn't contain any fasta file"
                return 1
            fi
        else
            rm -f "$genome_zip"

            if [ -d "$genome_tmp_unzipped_folder" ]; then
                rm -rf "$genome_tmp_unzipped_folder"
            fi

            echo "Failed to unzip the downloaded dataset"
            return 1
        fi
    else
        if [ -f "$genome_zip" ]; then
            rm -f "$genome_zip"
        fi

        echo "Failed to download genome using datasets command"
        return 1
    fi

    # Extract only relevant portion of fasta file based on description (case-insensitive)
    local best_match_header
    best_match_header=$(grep -i "^>" "$fasta_file" | grep -i "$genome_description" | head -n 1)
    
    if [ -z "$best_match_header" ]; then
        for word in $genome_description; do
            if [ ${#word} -gt 2 ]; then
                best_match_header=$(grep -i "^>" "$fasta_file" | grep -i "$word" | head -n 1)
                if [ -n "$best_match_header" ]; then
                    break
                fi
            fi
        done
    fi
    
    if [ -z "$best_match_header" ]; then
        best_match_header=$(grep "^>" "$fasta_file" | head -n 1)
        # echo "Warning: No matching sequence found for description '$genome_description', using first sequence"
        # echo "All sequence headers in $fasta_file:"
        # grep "^>" "$fasta_file"
        # echo "Comparing genome_description: '$genome_description'"
    fi
    
    if [ -n "$best_match_header" ]; then
        local sequence_id
        sequence_id=$(echo "$best_match_header" | sed 's/^>//' | awk '{print $1}')
        
        # Use sed to extract from matching header to next header (or end of file)
        sed -n "/^>$sequence_id/,/^>/p" "$fasta_file" | sed '$d' > "$genome_file"
        
        # If the last sequence in file, sed removes the last line incorrectly, so handle that case
        if [ ! -s "$genome_file" ]; then
            sed -n "/^>$sequence_id/,\$p" "$fasta_file" > "$genome_file"
        fi
        
        if [ -s "$genome_file" ] && grep -q "^>" "$genome_file"; then
            rm -rf "$genome_tmp_unzipped_folder"

            return 0
        else
            rm -rf "$genome_tmp_unzipped_folder"

            echo "Failed to extract matching sequence from fasta file"
            return 1
        fi
    else
        rm -rf "$genome_tmp_unzipped_folder"

        echo "No sequence headers found in fasta file"
        return 1
    fi
}

extract_genbank_from_csv() {
    local id="$1"

    local row
    if row=$(find_row_of_id "$id"); then
        echo "$row" | awk -F',' '{print $4}' | xargs
        return 0
    else
        return 1
    fi
}

extract_description_from_csv() {
    local id="$1"

    local row
    if row=$(find_row_of_id "$id"); then
        echo "$row" | awk -F',' '{print $3}' | xargs
        return 0
    else
        return 1
    fi
}

find_row_of_id() {
    local id="$1"
    
    local row
    row=$(awk -F',' -v target_id="$id" '$1 == target_id {print; found=1} END {if(!found) exit 1}' "$GENOME_CRISPR_CSV")
    
    if [ $? -eq 0 ] && [ -n "$row" ]; then
        echo "$row"
        return 0
    else
        echo "ID '$id' not found in CSV file"
        return 1
    fi
}

extract_crispr_sequences() {
    local success_count=0
    local total_count=0
    local successful_ids=""
    local pids=()
    local id_array=($BUCKET_IDS)

    for id in $BUCKET_IDS; do
        if [ -n "$id" ]; then
            ((total_count++))
            (
                if extract_single_crispr_sequences "$id" 2>/dev/null; then
                    exit 0
                else
                    exit 1
                fi
            ) &
            pids+=($!)
        fi
    done

    local i=0
    for pid in "${pids[@]}"; do
        wait "$pid"
        if [ $? -eq 0 ]; then
            ((success_count++))
            successful_ids="$successful_ids ${id_array[$i]}"
        fi
        ((i++))
    done

    BUCKET_IDS=$(echo "$successful_ids" | xargs)

    if [ "$success_count" -ne 0 ]; then
        echo "    ▸ ✔️  Extracted ${success_count}/${total_count} CRISPR sequences"
        return 0
    else
        echo "    ▸ ❌  Extracted ${success_count}/${total_count} CRISPR sequences"
        return 1
    fi
}

extract_single_crispr_sequences() {
    local id="$1"
    local genome_file="${GENOMES_FOLDER}/${id}.fasta"

    if [ -f "$genome_file" ]; then
        local genome_sequence
        if genome_sequence=$(grep -v ">" "$genome_file" | tr -d '\n'); then
            local crispr_regions
            if crispr_regions=$(extract_crispr_region_from_csv "$id"); then
                local output_file="${EXPECTED_CRISPRS_FOLDER}/${id}.txt"

                > "$output_file"

                # crispr_regions format: (1524978_148);(2004197_104);...
                echo "$crispr_regions" | sed -E 's/\(([0-9]+)_([0-9]+)\)/\1 \2/g' | tr ';' '\n' | while read start length; do
                    local end=$((start + length - 1))

                    local sequence
                    if sequence=$(echo "$genome_sequence" | cut -c "${start}-${end}"); then
                        echo "$sequence" >> "$output_file"
                    else
                        echo "Invalid start and end position to extract crispr sequence"
                        return 1
                    fi
                done

                return 0
            else
                echo "Unable to extract crispr region from csv file"
                return 1
            fi
        else
            echo "Unable to extract genome content from file"
            return 1
        fi
    else
        echo "Unable to find matching genome fasta file"
        return 1
    fi
}

extract_crispr_region_from_csv() {
    local id="$1"

    local row
    if row=$(find_row_of_id "$id"); then
        echo "$row" | awk -F',' '{print $6}' | xargs
        return 0
    else
        return 1
    fi
}

reduce_genomes_to_crispr_region() {
    local success_count=0
    local total_count=0
    local successful_ids=""
    local pids=()
    local id_array=($BUCKET_IDS)

    for id in $BUCKET_IDS; do
        if [ -n "$id" ]; then
            ((total_count++))
            (
                if reduce_single_genome_to_crispr_region "$id" 2>/dev/null; then
                    exit 0
                else
                    exit 1
                fi
            ) &
            pids+=($!)
        fi
    done

    local i=0
    for pid in "${pids[@]}"; do
        wait "$pid"
        if [ $? -eq 0 ]; then
            ((success_count++))
            successful_ids="$successful_ids ${id_array[$i]}"
        fi
        ((i++))
    done

    BUCKET_IDS=$(echo "$successful_ids" | xargs)

    if [ "$success_count" -ne 0 ]; then
        echo "    ▸ ✔️  Reduced ${success_count}/${total_count} genomes"
        return 0
    else
        echo "    ▸ ❌  Reduced ${success_count}/${total_count} genomes"
        return 1
    fi
}

reduce_single_genome_to_crispr_region() {
    local id="$1"
    local genome_file="${GENOMES_FOLDER}/${id}.fasta"
    local padding_to_include_surrounding_crispr_region=10000

    local full_genome_file_backup="${GENOMES_FOLDER}/${id}_full_genome.fasta"
    cp "$genome_file" "$full_genome_file_backup" >/dev/null 2>&1

    if [ ! -f "$genome_file" ]; then
        echo "Error: Genome file '$genome_file' does not exist for ID: $id"
        return 1
    fi

    local crispr_regions
    crispr_regions=$(extract_crispr_region_from_csv "$id")
    if [ $? -ne 0 ]; then
        echo "Unable to find the crispr region property for id $id"
        return 1
    fi

    # We assume the genome file has only one ">CP123 ..." description
    local header=$(head -n 1 "$genome_file")
    local genome_sequence=$(grep -v ">" "$genome_file" | tr -d '\n')
    local base_pairs_per_line=$(echo "$genome_sequence" | fold -w 80 | sed -n '2p' | wc -c)
    local genome_length=${#genome_sequence}
    
    local reduced_sequence=""
    while read start length; do
        if [ -n "$start" ] && [ -n "$length" ]; then
            local end=$((start + length - 1))
            local padded_start=$((start - $padding_to_include_surrounding_crispr_region))
            local padded_end=$((end + $padding_to_include_surrounding_crispr_region))
            
            # Boundaries
            [ $padded_start -lt 1 ] && padded_start=1
            [ $padded_end -gt $genome_length ] && padded_end=$genome_length
            
            local region_sequence=$(echo "$genome_sequence" | cut -c "${padded_start}-${padded_end}")
            reduced_sequence="${reduced_sequence}${region_sequence}"
        fi
    done < <(echo "$crispr_regions" | sed -E 's/\(([0-9]+)_([0-9]+)\)/\1 \2/g' | tr ';' '\n')
    
    echo "$header" > "$genome_file"
    echo "$reduced_sequence" | fold -w "$base_pairs_per_line" >> "$genome_file"
    
    return 0
}

generate_artificial_reads() {
    local success_count=0
    local total_count=0
    local successful_ids=""
    local pids=()
    local id_array=($BUCKET_IDS)

    for id in $BUCKET_IDS; do
        if [ -n "$id" ]; then
            ((total_count++))
            (
                if generate_artificial_reads_for_id "$id" 2>/dev/null; then
                    exit 0
                else
                    exit 1
                fi
            ) &
            pids+=($!)
        fi
    done

    local i=0
    for pid in "${pids[@]}"; do
        wait "$pid"
        if [ $? -eq 0 ]; then
            ((success_count++))
            successful_ids="$successful_ids ${id_array[$i]}"
        fi
        ((i++))
    done

    BUCKET_IDS=$(echo "$successful_ids" | xargs)

    if [ "$success_count" -ne 0 ]; then
        echo "    ▸ ✔️  Generated ${success_count}/${total_count} reads"
        return 0
    else
        echo "    ▸ ❌  Generated ${success_count}/${total_count} reads"
        return 1
    fi
}

generate_artificial_reads_for_id() {
    local id="$1"
    local genome_file="${GENOMES_FOLDER}/${id}.fasta"
    local reads_prefix="${READS_FOLDER}/${id}_reads"

    local coverage=20
    local read_length=125 # HiSeq 125, MiSeq 300, NextSeq 300, NovaSeq 150

    local n_reads
    if [ "$QUICK" = true ]; then
        local file_size
        file_size=$(stat -c%s "$genome_file")

        n_reads=$((file_size * $coverage / $read_length))
    else
        local genome_size
        genome_size=$(extract_genome_size "$id")

        n_reads=$((genome_size * $coverage / $read_length))
        if [ $? -ne 0 ]; then
            echo "Failed to extract the genome size from the CSV"
            return 1
        fi
    fi
    if [ "$n_reads" -le 0 ]; then
        echo "Failed trying to generate $n_reads reads"
        return 1
    fi

    if iss generate --genomes "$genome_file" --model hiseq --n_reads $n_reads --output "$reads_prefix" --cpus 1 >/dev/null 2>&1; then
        local vcf_pattern="${reads_prefix}*.vcf"
        if compgen -G "$vcf_pattern" > /dev/null 2>&1; then
            rm -f $vcf_pattern
        fi

        txt_pattern="${reads_prefix}*.txt"
        if compgen -G "$txt_pattern" > /dev/null 2>&1; then
            rm -f $txt_pattern
        fi

        return 0
    else
        echo "Failed to execute iss generate command"
        return 1
    fi
}

extract_genome_size() {
    local id="$1"

    local row
    if row=$(find_row_of_id "$id"); then
        echo "$row" | awk -F',' '{print $2}' | xargs
        return 0
    else
        return 1
    fi
}

run_mcaat_benchmark() {
    local success_count=0
    local total_count=0
    local successful_ids=""

    for id in $BUCKET_IDS; do
        if [ -n "$id" ]; then
            ((total_count++))

            if run_single_mcaat_benchmark "$id"; then
                ((success_count++))
                successful_ids="$successful_ids $id"
            fi
        fi
    done

    BUCKET_IDS=$(echo "$successful_ids" | xargs)

    if [ "$success_count" -ne 0 ]; then
        echo "    ▸ ✔️  Generated ${success_count}/${total_count} reads"
        return 0
    else
        echo "    ▸ ❌  Generated ${success_count}/${total_count} reads"
        return 1
    fi
}

run_single_mcaat_benchmark() {
    local id="$1"
    local r1_file="${READS_FOLDER}/${id}_reads_R1.fastq"
    local r2_file="${READS_FOLDER}/${id}_reads_R2.fastq"
    local crispr_sequences_file="${EXPECTED_CRISPRS_FOLDER}/${id}.txt"
    local result_file="${RESULTS_FOLDER}/${id}.txt"

    if "$MCAAT_BINARY" --input-files "$r1_file" "$r2_file" --ram "$RAM" --threads "$THREADS" --benchmark "$crispr_sequences_file" > "$result_file" 2>&1; then
        return 0
    else
        echo "Failed to execute mcaat binary with arguments"
        return 1
    fi
}

main() {
    echo "----------------------------------------------------"
    echo "🔹SETUP🔹"
    echo "----------------------------------------------------"

    verify_requirements
    verify_arguments "$@"
    setting_up_environment

    echo ""
    echo "----------------------------------------------------"
    echo "🔹RUN🔹"
    echo "----------------------------------------------------"
    
    local samples_done=0
    local max_bucket_size=1
    if [ "$QUICK" = true ]; then
        max_bucket_size=$THREADS
    fi

    while [ $samples_done -lt $NUM_SAMPLES ]; do
        local bucket_size=$((max_bucket_size > NUM_SAMPLES - samples_done ? NUM_SAMPLES - samples_done : max_bucket_size))
        BUCKET_IDS=$(echo "$IDS" | cut -d' ' -f$((samples_done + 1))-$((samples_done + bucket_size)))

        if [ "$bucket_size" -eq 1 ]; then
            echo "  ▸ Running id: $BUCKET_IDS"
        else
            echo "  ▸ Running Bucket $((samples_done + bucket_size))/$NUM_SAMPLES, using ids: $BUCKET_IDS"
        fi

        download_genomes || continue
        extract_crispr_sequences || continue
        
        if [ "$QUICK" = true ]; then
            reduce_genomes_to_crispr_region || continue
        fi
        
        generate_artificial_reads || continue
        run_mcaat_benchmark
        
        samples_done=$((samples_done + bucket_size))
    done
}

main "$@"
