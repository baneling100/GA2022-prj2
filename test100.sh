for i in {1..100}; do
	./run.sh
	cat maxcut.out >> results.out
done
