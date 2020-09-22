# Start container with volume attached
docker run --rm -d -it --name iris3-workflow-env -v c:/Users/flyku/Documents/GitHub/iris3-backend/tmp:/iris3 wangcankun100/iris3-workflow-env bin/bash
docker run --rm -d --name iris3-workflow-env -v /var/www/html/iris3/test:/iris3 wangcankun100/iris3-workflow-env

# Get container ID
docker ps -aqf "name=iris3-workflow-env"

# Enter docker container
docker exec -it iris3-workflow-env bash

# Stop container
docker stop iris3-workflow-env

# Start redis with ./tmp as persistent storage

docker run --name iris3-redis -v c:/Users/flyku/Documents/GitHub/iris3-backend/tmp:/data -d redis redis-server --appendonly yes