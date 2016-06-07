#
# Copyright (C) Telecom ParisTech
# 
# This file must be used under the terms of the CeCILL. This source
# file is licensed as described in the file COPYING, which you should
# have received as part of this distribution. The terms are also
# available at:
# http://www.cecill.info/licences/Licence_CeCILL_V1.1-US.txt
#

SHELL		:= /bin/bash

CC		:= gcc
CFLAGS		:= -Wall -c -O3
INCLUDES	:= -I./include
LD		:= gcc
LDFLAGS		:=
LIBS		:= -lm

OBJS		:= $(patsubst %.c,%.o,pa.c $(wildcard src/*.c))
EXECS		:= pa
GNUPLOT_DATA	:= dpa.dat average.dat
GNUPLOT_CMD	:= dpa.cmd average.cmd

.PHONY: help all clean ultraclean check archive

define HELP_message
Type:
  <make> or <make help> to get this help message
  <make all> to build everything
  <make check> to check the compliance of your code with specifications
  <make archive> to create the archive of your work
  <make clean> to clean a bit
  <make ultraclean> to really clean
endef
export HELP_message

help:
	@echo "$$HELP_message"

all: $(EXECS)

pa: $(OBJS)
	$(LD) $(LDFLAGS) $^ -o $@ $(LIBS)

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) $< -o $@

define NOREPORT_message
Report not found. Please write your report, name it report.txt (plain text
report) or report.pdf (PDF format) and run the script again. Note: if both
formats are found, only the PDF version will be used. Exiting..."
endef
export NOREPORT_message

define NOSTUDENT1_message
Student #1 cannot be undefined. Please edit the team.txt file and specify
at least student1. Exiting...
endef
export NOSTUDENT1_message

define CHECK_message
################################################################################
Please check the following information and, if they are correct, send an
un-signed and un-encrypted e-mail to:
  renaud.pacalet@telecom-paristech.fr
with subject exactly:
  HWSec-PA
and with the generated HWSec-PA.tgz archive attached. The e-mail body and other
attachments will be ignored.
################################################################################
endef
export CHECK_message

archive:
	@MandatoryFilesList="team.txt pa.c" && \
	for f in $$MandatoryFilesList; do \
		if [ ! -f $$f ]; then \
			echo "$$f file not found. Exiting..." && \
			exit -1; \
		fi; \
	done && \
	if [ -f report.pdf ]; then \
		report=report.pdf; \
	elif [ -f report.txt ]; then \
		report=report.txt; \
	else \
		echo "$$NOREPORT_message" && \
		exit -1; \
	fi && \
	. ./team.txt && \
	if [ -z "$$student1" ]; then \
		echo "$$NOSTUDENT1_message" && \
		exit -1; \
	fi && \
	if [ -f HWSec-PA.tgz ]; then \
		echo "HWSec-PA.tgz archive already exists. Rename or delete it first. Exiting..." && \
		exit -1; \
	fi && \
	tar hzcf HWSec-PA.tgz "$$report" $$MandatoryFilesList && \
	echo "$$CHECK_message" && \
	echo "Student #1: $$student1" && \
	if [ -z "$$email1" ]; then \
		echo "Email #1:   none (results will be e-mailed to sender)"; \
	else \
		echo "Email #1:   $$email1"; \
	fi && \
	if [ ! -z "$$student2" ]; then \
		echo "Student #2: $$student2" && \
		if [ -z "$$email2" ]; then \
			echo "Email #2:   none"; \
		else \
			echo "Email #2:   $$email2"; \
		fi; \
	fi && \
	echo "Report:     $$report"

check: pa pa.dat
	@d=`mktemp --directory --tmpdir=.` && \
	cd $$d && \
	k=`../pa ../pa.dat 10 2> /dev/null` && \
	cd .. && \
	rm -rf $$d && \
	if [[ ! "$$k" =~ ^0x[0-9a-fA-F]{12}$$ ]]; then \
		echo "pa.c invalid output:" && \
		echo "$$k" && \
		exit 1; \
	fi && \
	echo && \
	echo "Check OK. You can now make the archive."

clean:
	rm -f $(OBJS) $(GNUPLOT_DATA) $(GNUPLOT_CMD)

ultraclean: clean
	rm -rf $(EXECS)
