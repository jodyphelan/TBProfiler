from .text import write_text
from .pdf import write_pdf
from datetime import datetime
from pathogenprofiler import infolog
import json

def write_outputs(args,results):
    infolog("\nWriting outputs")
    infolog("---------------")
    json_output = args.dir+"/results/"+args.prefix+".results.json"
    text_output = args.dir+"/results/"+args.prefix+".results.txt"
    csv_output = args.dir+"/results/"+args.prefix+".results.csv"
    pdf_output = args.dir+"/results/"+args.prefix+".results.pdf"
    extra_columns = [x.lower() for x in args.add_columns.split(",")] if args.add_columns else []
    results["timestamp"] = datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    infolog(f"Writing json file: {json_output}")
    json.dump(results,open(json_output,"w"))
    extra_columns = [x.lower() for x in args.add_columns.split(",")] if args.add_columns else []
    if "pdf" in vars(args) and args.pdf:
        infolog(f"Writing pdf file: {pdf_output}")
        write_pdf(results,args.conf,pdf_output)
    if args.txt:
        infolog(f"Writing text file: {text_output}")
        write_text(results,args.conf,text_output,extra_columns,reporting_af=args.reporting_af,sep="\t")
    if args.csv:
        infolog(f"Writing csv file: {csv_output}")
        write_text(results,args.conf,csv_output,extra_columns,reporting_af=args.reporting_af,sep=",")
