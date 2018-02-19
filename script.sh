#!/bin/bash

#run_time=5
mkdir -p sana_results_sqrt


#mkdir -p "$i"_run
for run_time in 10
do
	mkdir -p "$run_time"min

	for base_network in SPombe yeast
	do
		for network in M0.05 M0.15 M0.25 M0.05Rewire0.01 M0.15Rewire0.01 M0.25Rewire0.01
		do
			
			for i in {1..1}
			do				
				result_path="$PWD"/"$base_network"_"$base_network$network"_"$i"
				mkdir -p "$result_path"
				
				for method in sana sana_improved
				do
					./sana -method $method -g1 "$base_network$network" -g2 "$base_network" -t $run_time  
					mv -v "$PWD"/sana.out "$result_path"/"$method".out
					mv -v "$PWD"/sana.out.align "$result_path"/"$method".out.align
					for csv in "$PWD"/*"$base_network".csv
					do
						mv ${csv} "$result_path"
					done
				done 
				
				
			done
			
		done
		for result_file in "$PWD"/"$base_network"_*
		do
			mv $result_file "$run_time"min
		done
	done

	for base_network in SCerevisiae
	do
		for network in M0.25Rewire0.01
		do
			
			for i in {1..1}
			do				
				result_path="$PWD"/"$base_network"_"$base_network$network"_"$i"
				mkdir -p "$result_path"
				
				for method in sana sana_improved
				do
					./sana -method $method -g1 "$base_network$network" -g2 "$base_network" -t $run_time  
					mv -v "$PWD"/sana.out "$result_path"/"$method".out
					mv -v "$PWD"/sana.out.align "$result_path"/"$method".out.align
					for csv in "$PWD"/*"$base_network".csv
					do
						mv ${csv} "$result_path"
					done
				done 
				
				
			done
			
		done
		for result_file in "$PWD"/"$base_network"_*
		do
			mv $result_file "$run_time"min
		done
	done

	for base_network in DMelanogaster
	do
		for network in M0.25Rewire0.01
		do
			
			for i in {1..2}
			do				
				result_path="$PWD"/"$base_network"_"$base_network$network"_"$i"
				mkdir -p "$result_path"
				
				for method in sana sana_improved
				do
					./sana -method $method -g1 "$base_network$network" -g2 "$base_network" -t $run_time  
					mv -v "$PWD"/sana.out "$result_path"/"$method".out
					mv -v "$PWD"/sana.out.align "$result_path"/"$method".out.align
					for csv in "$PWD"/*"$base_network".csv
					do
						mv ${csv} "$result_path"
					done
				done 
				
				
			done
			
		done
		for result_file in "$PWD"/"$base_network"_*
		do
			mv $result_file "$run_time"min
		done
	done

	#mv "$run_time"min "$i"_run
done

for run in "$PWD"/*_min
do		
	mv $run sana_results_sqrt
done

 
