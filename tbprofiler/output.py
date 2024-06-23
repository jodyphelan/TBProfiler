from .text import write_text
from .docx import write_docx
import logging
from .models import ProfileResult

def write_outputs(args,result: ProfileResult,template_file: str = None):
    logging.info("Writing outputs")

    json_output = args.dir+"/results/"+args.prefix+".results.json"
    text_output = args.dir+"/results/"+args.prefix+".results.txt"
    csv_output = args.dir+"/results/"+args.prefix+".results.csv"
    docx_output = args.dir+"/results/"+args.prefix+".results.docx"

    logging.info(f"Writing json file: {json_output}")
    open(json_output,"w").write(result.model_dump_json(indent=4))

    if args.txt:
        logging.info(f"Writing text file: {text_output}")
        write_text(result,args.conf,text_output,sep="\t",template_file=template_file)
    if args.csv:
        logging.info(f"Writing csv file: {csv_output}")
        write_text(result,args.conf,csv_output,sep=",",template_file = template_file)
    if args.docx:
        logging.info(f"Writing docx file: {docx_output}")
        if args.docx_plugin:
            args.docx_plugin = args.plugins[args.docx_plugin]
        write_docx(result,args.conf,docx_output,template_file = args.docx_template, plugin = args.docx_plugin)