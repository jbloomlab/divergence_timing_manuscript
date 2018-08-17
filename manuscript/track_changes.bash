#!/bin/bash
OLD="manuscript_original_submission"
NEW="manuscript"
TRACK="diff"

latexdiff "$OLD.tex" "$NEW.tex" > "$TRACK.tex"

# color additions green rather than blue
sed -i -e 's/{blue}/{green}/g' "$TRACK.tex"

pdflatex $TRACK && bibtex $TRACK && pdflatex $TRACK && pdflatex $TRACK
