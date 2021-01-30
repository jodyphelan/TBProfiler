

css_text = """
#logo {
  height: 20px;
  width: 100px;
}

article {
    height: 297mm;
}
body {
        /* font-family: -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,"Helvetica Neue",Arial,sans-serif,"Apple Color Emoji","Segoe UI Emoji","Segoe UI Symbol","Noto Color Emoji"; */
    font-size: 10px;
    font-weight: 400;
    line-height: 1.5;
    color: #212529;
    text-align: left;
}

header {
    margin: 10px;
}

.card {
    margin: 10px;
}

.card-header {
    background-color: rgba(0,0,0,.03);
    padding:.75rem 1.25rem;
    margin-bottom:0;
    border-bottom:1px solid rgba(0,0,0,.125);
}

.card-body {
    flex: 1 1 auto;
    padding:1.25rem;
}

.table {
    margin-bottom: 1rem;
    background-color: transparent;
}

.table td {
    padding: .75rem;
    vertical-align: top;
    border-top: 1px solid #dee2e6;
}

table {
    border-collapse: collapse;
}

thead {
    border-bottom: 2px solid #dee2e6;
}
"""


html_text = """
<html>
<head>
    <link rel="stylesheet" href="style.css">
</head>
<article class="">
    <header>
        <div class="">
          <img id="logo" src="https://raw.githubusercontent.com/jodyphelan/jodyphelan.github.io/master/img/tb-profiler-logo-rectangle.png" alt="error"></img>
        </div>
    </header>

	<div class="card border-dark mb-3">
		<div class="card-header"><b>Summary</b></div>
		<div class="card-body text-dark">
			<p><b>Sample name:</b> {{ result["id"] }}</p>
			<p><b>Analysis timestamp:</b> {{ result["timestamp"] }}</p>
			<p><b>Strain:</b> {{ result["sublin"] }}</p>
			<p><b>Drug-resistance:</b> {{result["drtype"]}}</p>
            <p><b>Comments:</b></p>
            <div style="text-align:right; margin-right: 100px;">
                <p>Signed:</p>
                <p>Date:</p>
		</div>
	</div>


    <div class="card border-dark">
        <div class="card-header"><b>Drug resistance Mutations:</b> This table reports
            mutations found in candidate resistance genes which have been
            associated with drug resistance
        </div>
        <div class="card-body">
            {% if result["qc"]["gene_coverage"]|length==0 %}
                <div class="">
                    All genes have sufficient coverage
                </div>
            {% else %}
                <table class="table">
                    <thead>
                        <tr>
                            <th scope="col">Gene</th>
                            <th scope="col">Mutation</th>
                            <th scope="col">Type</th>
                            <th scope="col">Drugs</th>
                            <th scope="col">Estimated frequency (%)</th>
                        </tr>
                    </thead>
                    <tbody>
                    {% for var in result["dr_variants"] %}
                        {% set percent = var["freq"] * 100 %}
                        <tr>
                            <td>{{ var["gene"] }}</td>
                            <td>{{ var["change"] }}</td>
                            <td>{{ var["type"] }}</td>
                            <td>{{ var["drugs"] }}</td>
                            <td>{{ percent|int }}</td>
                        </tr>
                    {% endfor %}
                    </tbody>
                </table>
            {% endif %}
        </div>
    </div>


</article>

<article class="">
    <header>
        <div class="">
          <img id="logo" src="https://raw.githubusercontent.com/jodyphelan/jodyphelan.github.io/master/img/tb-profiler-logo-rectangle.png" alt="error"></img>
        </div>
    </header>

    <div class="card border-dark">
        <div class="card-header"><b>Other Mutations:</b> This table reports
            non-synonymous mutations found in candidate resistance genes which have not been
            associated with drug resistance</div>
        <div class="card-body">
            {% if result["qc"]["gene_coverage"]|length==0 %}
                <div class="">
                    All genes have sufficient coverage
                </div>
            {% else %}
                <table class="table">
                    <thead>
                        <tr>
                            <th scope="col">Gene</th>
                            <th scope="col">Mutation</th>
                            <th scope="col">Type</th>
                            <th scope="col">Estimated frequency (%)</th>
                        </tr>
                    </thead>
                    <tbody>
                    {% for var in result["other_variants"] %}
                        {% set percent = var["freq"] * 100 %}
                        <tr>
                            <td>{{ var["gene"] }}</td>
                            <td>{{ var["change"] }}</td>
                            <td>{{ var["type"] }}</td>
                            <td>{{ percent|int }}</td>
                        </tr>
                    {% endfor %}
                    </tbody>
                </table>
            {% endif %}
        </div>
    </div>

</article>


<article>
    <header>
        <div class="">
          <img id="logo" src="https://raw.githubusercontent.com/jodyphelan/jodyphelan.github.io/master/img/tb-profiler-logo-rectangle.png" alt="error"></img>
        </div>
    </header>
    <div class="card border-dark">
        <div class="card-header"><b>Low coverage genes</b></div>
        <div class="card-body">
            {% if result["qc"]["gene_coverage"]|length==0 %}
                <div class="">
                    All genes have sufficient coverage
                </div>
            {% else %}
                <table class="table">
                    <thead>
                        <tr>
                            <th scope="col">Gene</th>
                            <th scope="col">Percent </th>
                            <th scope="col">Cutoff</th>
                        </tr>
                    </thead>
                    <tbody>
                    {% for gene in result["qc"]["gene_coverage"] %}
                        {% set percent = gene["fraction"] * 100 %}
                        <tr>
                            <td>{{ gene["gene"] }}</td>
                            <td>{{ percent|int }}</td>
                            <td>{{ gene["cutoff"] }}</td>
                        </tr>
                    {% endfor %}
                    </tbody>
                </table>
            {% endif %}
        </div>
    </div>

    <div class="card border-dark mb-3 shadow">
		<div class="card-header bg-dark text-white text-center"><b>Pipeline</b></div>
		<div class="card-body text-dark">
			<p><b>TB-Profiler version:</b> {{ result["tbprofiler_version"] }}</p>
            <p><b>Database name:</b> {{ result["db_version"]["name"] }}</p>
            <p><b>Database version:</b> {{ result["db_version"]["commit"] }}</p>
            <p><b>Mapping software:</b> {{ result["pipeline"]["mapper"] }}</p>
            <p><b>Variant calling software:</b> {{ result["pipeline"]["variant_caller"] }}</p>
		</div>
	</div>
</article>


</html>
"""



def write_pdf(results,conf,outfile):
    from jinja2 import Environment, FileSystemLoader
    from weasyprint import HTML, CSS

    env = Environment(loader=FileSystemLoader('.'))
    template = env.from_string(html_text)

    gene_cov = []
    for gene in results["qc"]["gene_coverage"]:
        if gene["fraction"]==0:
            continue
        gene_cov.append(gene)

    dr_variants = []
    for var in results["dr_variants"]:
        var["drugs"] = ", ".join([d["drug"] for d in var["drugs"]])
        dr_variants.append(var)

    other_variants = []
    for var in results["other_variants"]:
        if "synonymous" in var["type"] or "stop_retained" in var["type"]:
            continue
        other_variants.append(var)

    results["dr_variants"] = dr_variants
    results["other_variants"] = other_variants
    results["qc"]["gene_coverage"] = gene_cov
    template_vars = {"result" : results}
    html_out = template.render(template_vars)
    HTML(string=html_out).write_pdf(outfile,stylesheets=[CSS(string=css_text)])
