testpypi_install: ## Install from testpypi
	virtualenv -p python3 ../testpypi
	../testpypi/bin/pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple masstodon

testpypi_upload: ## Upload repo to testpypi
	twine upload -r test dist

pypi_upload: ## upload repo to pypi
	twine upload -r pypi --skip-existing dist/*

pypi_install: ## Install from testpypi
	virtualenv -p python3 ../pypi
	../pypi/bin/pip install masstodon



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
