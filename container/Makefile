master-kill:
	docker kill m1
master-rm:
	docker rm m1
master-start:
	docker run --net rclust --ip 172.18.0.10 --name m1 -it gveeden/rcluster:Master
master-build:
	docker build -t gveeden/rcluster:Master /home/gerald/Dev/rParallel/Master/

workers-kill:
	docker kill w1 && docker kill w2 && docker kill w3
workers-rm:
	docker rm w1 && docker rm w2 && docker rm w3
workers-start:
	docker run -d -P --net rclust --ip 172.18.0.2 --name w1 -it gveeden/rcluster:Workers && docker run -d -P --net rclust --ip 172.18.0.3 --name w2 -it gveeden/rcluster:Workers && docker run -d -P --net rclust --ip 172.18.0.4 --name w3 -it gveeden/rcluster:Workers
workers-build:
	docker build -t gveeden/rcluster:Workers /home/gerald/Dev/rParallel/Workers/
