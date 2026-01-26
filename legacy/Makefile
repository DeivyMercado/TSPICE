#VARIABLES
SHELL:=/bin/bash
BRANCH=$(shell bash bin/getbranch.sh)
VERSION=$(shell tail -n 1 .versions 2>/dev/null || echo "0.0.1")
COMMIT=[MAN] Maintenance
RELMODE=release
PYTHON=python3
PIP=pip3
PACKNAME=TSPICE

#Show
show:
	@echo "Version: $(VERSION)"
	@echo "GitHub branch: $(BRANCH)"

##Clean tasks
clean:cleancrap

cleanall:cleancrap cleanout cleandist

cleancrap:
	@echo "Cleaning crap..."
	@-find . -name "*~" -delete
	@-find . -name "#*#" -delete
	@-find . -name "__pycache__" -type d | xargs rm -fr
	@-find . -name ".ipynb_checkpoints" -type d | xargs rm -fr
	@-find . -name ".DS_Store" -delete

cleanout:
	@echo "Cleaning all compiled objects..."
	@-find . -name "*.pyc" -delete

cleandist:
	@-rm -rf dist/
	@-rm -rf build/
	@-rm -rf *.egg-info/

#Git task
addall:cleanall
	@echo "Adding..."
	@-git add -A .

commit:
	@echo "Committing..."
	@git commit -am "$(COMMIT)"
	@-git push origin $(BRANCH)

pull:
	@echo "Pulling new files..."
	@-git reset --hard HEAD
	@-git pull origin $(BRANCH)

#Package tasks
release:
	@echo "Releasing a new version..."
	@bash bin/release.sh $(RELMODE) $(VERSION)

install:
	@$(PIP) install -e .

test:
	@$(PYTHON) -c "import tspice; print('tspice imported successfully')"