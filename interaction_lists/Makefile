ruan_2015:
	@bash ruan_2015.commands

# After ChIA-PET pipelines
motifs:
	bash find_motifs.commands

thurman: thurman.final_interactions
thurman.final_interactions :
	bash ../scripts/extract_thurman_interactions.sh

fantom: fantom.final_interactions
fantom.final_interactions :
	bash ../scripts/extract_fantom_interactions.sh

# Must be done after motifs
combine_E_interactions: E.final_interactions
E.final_interactions :
	bash ../scripts/combine_E_interactions.sh

# NOTE: This includes ChIA-PET data only
# FANTOM and Thurman P2P interactions are separate
combine_P2P_interactions: P2P.final_interactions
P2P.final_interactions :
	bash ../scripts/combine_P2P_interactions.sh
