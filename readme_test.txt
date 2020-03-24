datafunk clean_names -i test_files/clean_test/caseData.csv -t country/region -o test_files/clean_test/filtered.csv
datafunk merge_fasta -f test_files/merge_test -i test_files/merge_test/metadata.csv -o test_files/merge_test/merged.fasta
datafunk remove_fasta -i test_files/filter_test/test_remove.fasta -f test_files/filter_test/test_filter.txt -o test_files/filter_test/filtered.fasta
