# Exemplo de Início Rápido

Execute o pipeline inteiro:

```
tb-profiler profile -1 /path/to/reads_1.fastq.gz -2 /path/to/reads_2.fastq.gz -p prefix
```

O argumento `-p` permite fornecer um prefixo para os arquivos de saída resultantes (output). Isso é útil quando você precisa executar mais de uma amostra. Isso armazenará BAM, VCF e arquivos de resultado nos respectivos diretórios. Os resultados são produzidos em formatos json e texto.

## Exemplo de execução

```
mkdir test_run; cd test_run
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR166/009/ERR1664619/ERR1664619_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR166/009/ERR1664619/ERR1664619_2.fastq.gz
tb-profiler profile -1 ERR1664619_1.fastq.gz -2 ERR1664619_2.fastq.gz -t 4 -p ERR1664619
cat results/ERR1664619.results.json
```

## Executando com um arquivo BAM existente

Usando a opção `-a`, você pode especificar o uso de um arquivo BAM existente em vez de arquivos fastq.

```
tb-profiler profile -a /path/to/bam -p test
```

!!! warning
    Atenção! Os arquivos BAM devem ter sido criados usando a versão do genoma como banco de dados que pode ser baixado aqui. Surpreendentemente, este genoma tem vários números de acesso (ASM19595v2, NC_000962.3, GCF_000195955.2, etc.). Se você acredita que sua referência é exatamente a mesma sequência (o comprimento deve ser 4411532), você pode criar um banco de dados com o mesmo nome de sequência usado em seu arquivo BAM. Por exemplo, se o nome da sequência for "NC_000962.3", você pode fazer isso executando o seguinte:

    ```
    tb-profiler update_tbdb --match_ref /path/to/your/reference.fasta
    ```