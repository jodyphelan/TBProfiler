# Gerando arquivos resumidos

Os resultados de várias execuções podem ser agrupados em uma tabela usando o seguinte comando:

```
tb-profiler collate
```

Isso criará automaticamente um número de arquivos de resultados agrupados de todos os arquivos de resultados individuais no diretório result. Se você gostaria de gerar este arquivo para um subconjunto das execuções, você pode fornecer uma lista com a mesma execução usando o sinalizador `--samples`. O prefixo para os arquivos de saída é tbprofiler por padrão, mas isso pode ser alterado com o sinalizador `--prefix`.

## Escrevendo seus próprios scripts de resumo

A função `collate` (agrupar) extrai as mutações de resistência às drogas e de linhagem; no entanto, você pode querer extrair mais recursos que estão presentes nos arquivos de resultados individuais no formato json. Criei um pequeno tutorial sobre como fazer isso [aqui](https://jodyphelan.gitbook.io/tb-profiler/writing-a-custom-collate-script)

## Arquivos do iTOL

Vários arquivos são produzidos pela função `tb-profile collate`. Entre eles estão vários arquivos de configuração que podem ser usados no iTOL para anotar árvores filogenéticas. Uma pequena árvore e arquivos de configuração foram adicionados no diretório example_data. Para usar, navegue até o site do iTOL e carregue o arquivo tbprofiler.tree usando o botão de upload na barra de navegação. Depois de fazer o upload, você será levado a uma visualização da árvore. Para adicionar a anotação, clique no botão '+' no canto inferior direito e selecione os arquivos de configuração do iTOL. Agora você deve ver uma figura semelhante à abaixo. As seguintes anotações estão incluídas:

* Linhagem
* Classes de resistência às drogas (sensível, droga resistente, MDR e XDR)
* A droga resistência por cada droga individualmente, onde os círculos preenchidos representam a resistência.

<img href="https://github.com/jodyphelan/TBProfiler/raw/docs/docs/docs/assets/images/itol_example.png">

