ruan_2015:
	@bash ruan_2015.commands

# After ChIA-PET pipelines
motifs:
	bash find_motifs.commands

thurman: thurman.interactions
thurman.interactions :
	bash ../scripts/extract_thurman_interactions.sh

fantom: fantom.interactions
fantom.interactions :
	bash ../scripts/extract_fantom_interactions.sh

# Must be done after motifs
combine_E_interactions: E.interactions
E.interactions :
	bash ../scripts/combine_E_interactions.sh

combine_P2P_interactions: P2P.interactions
P2P.interactions :
	bash ../scripts/combine_P2P_interactions.sh
