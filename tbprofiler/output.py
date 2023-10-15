from .text import write_text
# from .pdf import write_pdf
from .docx import write_docx
from datetime import datetime
import logging
import json
import tbprofiler as tbp

def write_outputs(args,results,template_file = None):
    logging.info("Writing outputs")

    json_output = args.dir+"/results/"+args.prefix+".results.json"
    text_output = args.dir+"/results/"+args.prefix+".results.txt"
    csv_output = args.dir+"/results/"+args.prefix+".results.csv"
    docx_output = args.dir+"/results/"+args.prefix+".results.docx"
    # tree_output = args.dir+"/results/"+args.prefix+".results.nwk"

    if "add_columns" not in vars(args):
        args.add_columns = None
    extra_columns = [x.lower() for x in args.add_columns.split(",")] if args.add_columns else []
    results["timestamp"] = datetime.now().strftime("%d-%m-%Y %H:%M:%S")

    logging.info(f"Writing json file: {json_output}")
    json.dump(results,open(json_output,"w"))

    if args.txt:
        logging.info(f"Writing text file: {text_output}")
        write_text(results,args.conf,text_output,extra_columns,sep="\t",template_file=template_file)
    if args.csv:
        logging.info(f"Writing csv file: {csv_output}")
        write_text(results,args.conf,csv_output,extra_columns,sep=",",template_file = template_file)
    if args.docx:
        logging.info(f"Writing docx file: {docx_output}")
        write_docx(results,args.conf,docx_output,template_file = args.docx_template)