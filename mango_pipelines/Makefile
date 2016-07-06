### Download mango
# Note that before it is run, you should make sure that the requisite packages are installed: https://github.com/dphansti/mango

install_mango:
	@git clone https://github.com/dphansti/mango.git
	@R CMD INSTALL --no-multiarch --with-keep.source mango

prepare_workspace:
	@bash setup/prepare_workspace.sh

### Generate the pipeline files for running the ENCODE chia-pet experiments through mango
ENCODE_pipelines:
	@bash setup/make_ENCODE_chiapet_pipelines.sh

### pipeline files should be run individually
