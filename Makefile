testpypi_upload: ## Upload repo to testpypi
	twine upload -r test dist

pypi_upload: ## Upload repo to pypi
	twine upload -r pypi --skip-existing dist/*

testpypi_install: ## Install from testpypi
	rm -rf ../testpypi || true
	virtualenv -p python3 ../testpypi
	../testpypi/bin/pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple masstodon

pypi_install: ## Install from testpypi
	rm -rf ../pypi || true
	virtualenv -p python3 ../pypi
	../pypi/bin/pip install masstodon

test_masstodon_cli: ## testing command line tools
	./bin/masstodon spectrum.txt -m 11=H-1C10Ag2 5N=H-1CAg 4C_carbo=HPO

# -----------------------------------------------------------
# -----  EVERYTHING BELOW THIS LINE IS NOT IMPORTANT --------
# -----       (Makefile helpers and decoration)      --------
# -----------------------------------------------------------
#
# Decorative Scripts - Do not edit below this line unless you know what it is

.PHONY: help
.DEFAULT_GOAL := help

NO_COLOR    = \033[0m
INDENT      = -30s
BIGINDENT   = -50s
GREEN       = \033[36m
BLUE        = \033[34m
DIM         = \033[2m
help:
	@printf '\n\n$(DIM)Commands:$(NO_COLOR)\n'
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "$(GREEN) % $(BIGINDENT)$(NO_COLOR) - %s\n", $$1, $$2}'
