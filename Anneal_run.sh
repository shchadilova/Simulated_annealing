for vr in 0.5
do 
for M in 5
do  
    mkdir ../Sim_anneal_data/Run2/$M/
    
    echo MF $vr $M
    ./Anneal_coherent $vr $M >>../Sim_anneal_data/Run2/$M/history_mf_$vr.dat
	cp alpha.dat ../Sim_anneal_data/Run2/$M/alpha_mf_$vr.dat
    echo Diag $vr $M
	./Anneal_gaussian $vr $M >>../Sim_anneal_data/Run2/$M/history_diag_$vr.dat
	cp alpha.dat ../Sim_anneal_data/Run2/$M/alpha_diag_$vr.dat
	cp gamma.dat ../Sim_anneal_data/Run2/$M/gamma_diag_$vr.dat
	
done 
done
