for i in {1..10}; do
	./run.sh
	cat maxcut.out >> results.out
done
