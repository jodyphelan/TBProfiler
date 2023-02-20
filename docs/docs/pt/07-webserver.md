# Servidor web

## LSHTM TB-Profiler webserver

Se você não tem acesso a um ambiente linux ou macOS, você ainda pode usar o tb-profiler usando nosso servidor web http://tbdr.lshtm.ac.uk/

Você pode carregar seus dados de sequenciamento de próxima geração (next generation sequencing) no formato fastQ. Você pode carregar um ou dois arquivos fastq (forward e reverse). Quando você faz upload de seus dados, a corrida recebe um ID exclusivo. Anote este ID, pois você precisará encontrar seus resultados posteriormente. O upload em lote de amostras também é possível.

## Configurando seu próprio servidor web

O código do servidor web está disponível [aqui](https://github.com/jodyphelan/tb-profiler-webserver). Isso ainda está em desenvolvimento inicial, mas você pode usá-lo para configurar sua própria instância do servidor. Para fazer isso, execute o seguinte código:

```
# Install libraries
git clone https://github.com/jodyphelan/tb-profiler-webserver.git
cd tb-profiler-webserver
python setup.py install

# Run flask
export FLASK_APP=tbprofiler_web
export FLASK_ENV=development
flask run

# Run rabbit-mq server
rabbitmq-server

# Run celery
celery -A tbprofiler_web.worker worker --loglevel=info --concurrency=1
```
