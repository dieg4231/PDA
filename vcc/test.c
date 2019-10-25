#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include  <signal.h>


void     INThandler(int);
 
FILE * fp;

int main(void) {
signal(SIGINT, INThandler);
    fp = popen("python simulator_node.py", "w");
    if (fp == NULL) {
        printf("popen error\n");
        exit(1);
    }

    float inc = 0;
    char buf[10];

    while(1) {
        fprintf(fp, "%f \n", inc);
        inc+=.1;
        fflush(fp);
	sleep(1);
    }


    return (0);
}



void  INThandler(int sig)
{
     signal(sig, SIG_IGN);
     printf("OUCH \n ");
     pclose(fp);
     exit(0);

}















