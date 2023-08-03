from pathogenprofiler import infolog
import importlib.util


def apply_rules(conf,results):
    spec = importlib.util.spec_from_file_location("tbprofiler.rules", conf['rules'])
    rules = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(rules)

    rule_names = [k for k in vars(rules) if k[0]!="_"]
    for r in rule_names:
        infolog(f"Applying rule: {r}")
        results = vars(rules)[r](results)

