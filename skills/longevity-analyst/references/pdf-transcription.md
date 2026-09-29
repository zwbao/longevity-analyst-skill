# Transcribing lab reports (PDF, photos)

The ground truth is what the page shows, read with multimodal vision page by page. A text layer or OCR may help you
find things. It never overrides what you see.

1. Read every page. For each lab row, write `marker,value,unit,ref_range,page` to a CSV, exactly as printed
   (`ref_range` is the printed reference interval, empty if none):
   the marker name as printed (Chinese and abbreviation both, e.g. `白蛋白 ALB`), the value with its sign
   (`<0.5` stays `<0.5`), the unit as printed. No conversion, no rounding.
2. A value you cannot read with confidence is left out and listed in `--reason`.
3. Numeric body measurements printed on the page (blood pressure, waist, weight, height) **are** rows too, because
   methods read them from the measurement table: `收缩压,136,mmHg,,1` and `舒张压,84,mmHg,1` for "血压 136/84 mmHg",
   `腰围,94,cm,1`. Yes/no facts (smoking, medication use, diabetes diagnosis) are not rows; record them with
   `la.py member <ws> key=value --source "<file> p<page>: <quote>"` when a method asks for them.
4. `la.py labs add <ws> --csv rows.csv --file <Fxxx> --reason "transcribed N rows from p1-2; unreadable: ..."`.
