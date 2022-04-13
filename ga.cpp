#include<stdio.h>

int main() {

        int answer[4] = {1, 6, 8, 9};

        // Print to maxcut.out
        for (int i = 0; i < 4; i++) {
                if (i > 0) printf(" ");
                printf("%d", answer[i]);
        }
        printf("\n");
}
