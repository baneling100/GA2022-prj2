for i in {1..500}; do
	make && ./ga < maxcut.in > maxcut.out
	cat maxcut.out >> results.out
done
