grep "n_pseudoaligned" mapping/*/run_info.json | awk '{print $1 $3}' | sed -e 's/mapping[/]//g' -e 's/[/]run_info.json:/\t/g' -e 's/,//g' > mapping.rate

