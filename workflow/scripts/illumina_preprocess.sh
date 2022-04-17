#!/usr/bin/env bash
set -xeo pipefail

DEFAULT_THREADS=1

function usage {
    echo "usage: illumina_preprocess.sh -r <STR> -i <FILE> -R <FILE> -o <FILE> [OPTIONS]"
    echo "   "
    echo "  -r                       : Run accession [REQUIRED]"
    echo "  -i                       : Run info TSV file [REQUIRED]"
    echo "  -R                       : Output report HTML file [REQUIRED]"
    echo "  -o                       : Output fastq file [REQUIRED]"
    echo "  -t                       : Number of fastp threads [default: $DEFAULT_THREADS]"
    echo "  -h | --help              : This message"
}

function parse_args {
    # positional args
    args=()

    # named args
    while [ "$1" != "" ]; do
        case "$1" in
            -r)
                run_acc="$2"
                shift
                ;;
            -i)
                run_info="$2"
                shift
                ;;
            -R)
                report="$2"
                shift
                ;;
            -o)
                output="$2"
                shift
                ;;
            -t)
                threads="$2"
                shift
                ;;
            -h | --help)
                usage
                exit
                ;;             # quit and show usage
            *) args+=("$1") ;; # if no match, add it to the positional args
        esac
        shift # move to next kv pair
    done

    # restore positional args
    #    set -- "${args[@]}"

    # set positionals to vars
    #  positional_1="${args[0]}"
    #  positional_2="${args[1]}"

    # validate required args
    if [[ -z "${run_acc}" || -z "${run_info}" || -z "${report}" || -z "${output}" ]]; then
        echo "Invalid arguments"
        usage
        exit
    fi

    # set defaults
    if [[ -z "$threads" ]]; then
        threads="$DEFAULT_THREADS"
    fi
}

function run {
    parse_args "$@"

    files_str=$(grep "$run_acc" "$run_info" | cut -f2)
    IFS=';' read -r -a files <<< "$files_str"
    n_files="${#files[@]}"

    input_arg=("-i ${files[0]}")
    if [ "$n_files" -eq 2 ]; then
        input_arg+=("-I ${files[1]}")
        args+=("--detect_adapter_for_pe")
    fi

    fastp -h "$report" -w $threads ${input_arg[*]} ${args[*]} | gzip -c > "$output"
}

run "$@"
