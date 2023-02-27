# Biblioteca das mutações

TB-Profiler vem com um banco de dados padrão. O desenvolvimento da biblioteca de mutação está hospedado no repositório tbdb. Visite este repositório se você deseja se envolver com o banco de dados ou modificar e criar o seu próprio.

## Por que existe um repositório github separado?

Com pipelines de análise praticamente padronizados, é evidente que a acurácia da predição é afetada principalmente pela biblioteca subjacente de mutações. À medida que novas evidências para a inclusão ou exclusão de mutações são geradas, há uma necessidade constante de atualizar e reavaliar a biblioteca de mutações. Além disso, é importante que o controle da biblioteca fique nas mãos dos usuários finais. Ao hospedar a biblioteca em um repositório separado (em vez de ocultado profundamente no código da ferramenta para criação de perfil), fica mais fácil descobrir exatamente quais mutações estão presentes. Além disso, o github tem uma série de recursos úteis que podem ser utilizados:

* Todas as alterações na biblioteca são rastreadas e podem ser facilmente revertidas para versões anteriores.
* Os usuários podem levantar questões ou discutir a biblioteca usando a seção de [problemas](https://github.com/jodyphelan/tbdb/issues) (issues) do repositório github.
* Diferentes versões da biblioteca podem ser mantidas usando [Forks](https://help.github.com/en/articles/fork-a-repo), permitindo que os usuários experimentem a biblioteca sem afetar o projeto principal. Essas alterações podem ser mescladas no repositório principal após as alterações serem revisadas.
* Vários usuários / desenvolvedores podem contribuir para a biblioteca.

!!! info
    tl;dr – Hospedando este separadamente torna mais fácil atualizar a biblioteca.

## Quer contribuir?

Se você acha que uma mutação deve ser removida ou adicionada, indique e questione [aqui](https://github.com/jodyphelan/tbdb/issues). Se você quiser ajudar na curadoria da biblioteca, deixe um comentário [aqui](https://github.com/jodyphelan/tbdb/issues/4)

## Adicionando/removendo mutações

As mutações podem ser adicionadas enviando uma pull request (PR) com o arquivo tbdb.csv modificado. Se a frase anterior não fazia sentido para você, você pode sugerir uma alteração usando um [problema](https://github.com/jodyphelan/tbdb/issues/4) e tentaremos ajudar. 