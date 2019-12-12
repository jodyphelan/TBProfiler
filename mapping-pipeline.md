# Mapping pipeline

First log into the server

```
ssh monica@10.18.0.21
```

Then go to the correct directory with the data

```text
cd /opt/storage5/monica/NGS_data
```



{% hint style="info" %}
If you have a new reference you have to index it using

```text
bwa index <reference_file>
```
{% endhint %}

 Perform mapping using the following code

```text
bwa mem -t 10 -R "@RG\tID:AGT12\tSM:ATG12\tPL:Illumina" Anopheles_gambiae.AgamP4.dna.toplevel.fa AGT12_1.fastq.gz AGT12_2.fastq.gz | samtools view -@ 10 -b - | samtools sort -@ 10 -o AGT12.bam -
```

Then index the bam file with

```text
samtools index AGT12.bam
```

