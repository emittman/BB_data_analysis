all: .makefile/datasets
all: .makefile/analyses
all: .makefile/data_3_10
all: .makefile/samples_extreme_relaxation_3_24

.makefile/datasets: .makefile/data_3_10

.makefile/analyses: .makefile/samples_extreme_relaxation_3_24

.makefile/data_3_10:
	${PARALLEL_REMAKE_RUNNER} Rscript -e 'parallelRemake:::process("$@", remake_file = "remake.yml", verbose = TRUE)'

.makefile/samples_extreme_relaxation_3_24: .makefile/data_3_10
.makefile/samples_extreme_relaxation_3_24:
	${PARALLEL_REMAKE_RUNNER} Rscript -e 'parallelRemake:::process("$@", remake_file = "remake.yml", verbose = TRUE)'

