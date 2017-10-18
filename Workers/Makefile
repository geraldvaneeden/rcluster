kill:
	docker kill w1 && docker kill w2 && docker kill w3
rm:
	docker rm w1 && docker rm w2 && docker rm w3
start:
	docker run -d -P --net rclust --ip 172.18.0.2 --name w1 -it gveeden/rcluster:Workers && docker run -d -P --net rclust --ip 172.18.0.3 --name w2 -it gveeden/rcluster:Workers && docker run -d -P --net rclust --ip 172.18.0.4 --name w3 -it gveeden/rcluster:Workers
build:
	docker build -t gveeden/rcluster:Workers /home/gerald/Dev/rParallel/Workers/
