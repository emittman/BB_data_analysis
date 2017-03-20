all: .makefile/datasets
all: .makefile/analyses
all: .makefile/data_2_21
all: .makefile/data_3_10
all: .makefile/samples_logodds_reduced_3_8_data_2_21
all: .makefile/samples_logodds_reduced_3_8_data_3_10

.makefile/datasets: .makefile/data_2_21
.makefile/datasets: .makefile/data_3_10

.makefile/analyses: .makefile/samples_logodds_reduced_3_8_data_2_21
.makefile/analyses: .makefile/samples_logodds_reduced_3_8_data_3_10

.makefile/data_2_21:
	${PARALLEL_REMAKE_RUNNER} Rscript -e 'parallelRemake:::process("$@", remake_file = "remake.yml", verbose = TRUE)'

.makefile/data_3_10:
	${PARALLEL_REMAKE_RUNNER} Rscript -e 'parallelRemake:::process("$@", remake_file = "remake.yml", verbose = TRUE)'

.makefile/samples_logodds_reduced_3_8_data_2_21: .makefile/data_2_21
.makefile/samples_logodds_reduced_3_8_data_2_21:
	${PARALLEL_REMAKE_RUNNER} Rscript -e 'parallelRemake:::process("$@", remake_file = "remake.yml", verbose = TRUE)'

.makefile/samples_logodds_reduced_3_8_data_3_10: .makefile/data_3_10
.makefile/samples_logodds_reduced_3_8_data_3_10:
	${PARALLEL_REMAKE_RUNNER} Rscript -e 'parallelRemake:::process("$@", remake_file = "remake.yml", verbose = TRUE)'

