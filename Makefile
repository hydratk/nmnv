.PHONY: clean

# Directories
REPORT_DIR = report
SECTIONS_DIR = $(REPORT_DIR)/sections

# File patterns
CLEAN_EXTENSIONS = *.aux *.fdb_latexmk *.fls *.log *.nav *.snm *.toc *.out *.pdf *.synctex.gz

clean:
	rm -f $(REPORT_DIR)/$(CLEAN_EXTENSIONS)
	rm -f $(SECTIONS_DIR)/*.aux
	rm -f $(REPORT_DIR)/project.synctex.gz
