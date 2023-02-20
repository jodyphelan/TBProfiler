# Atualizando

TB-Profiler está em constante desenvolvimento rápido. Se você planeja usar o programa em seu trabalho, certifique-se de usar a versão mais recente! Da mesma forma, o banco de dados não é estático e está sendo continuamente aprimorado, portanto, certifique-se de usar a versão mais recente. Se você usa TBProfiler em seu trabalho, indique a versão da ferramenta e do banco de dados, pois foram desenvolvidos de forma independente um do outro.

## Atualizando o banco de dados

Novas mutações / genes são adicionados periodicamente ao banco de dados. Execute o seguinte para se certificar de que está atualizado.

```
tb-profiler update_tbdb
```

## Refazendo o perfil rapidamente

Se você tem um novo banco de dados de mutação, mas nenhum novo gene foi adicionado, você pode rapidamente refazer o perfil de suas amostras executando o seguinte.

```
tb-profiler reprofile /path/to/result.json
```

Isso pode ser útil quando você mesmo adicionou algumas mutações ou tem certeza de que nenhum novo gene foi adicionado na atualização. Se você não tiver certeza, é mais seguro executar a etapa completa de criação de perfil novamente.

