#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <math.h>

#include <jack/jack.h>

// Include FFTW header
#include <complex.h> //needs to be included before fftw3.h for compatibility
#include <fftw3.h>

#define WINDOW_SIZE 1024
#define WINDOWS_PER_BUFF_0 6
#define WINDOWS_PER_BUFF_AB 4
#define NUM_CH  2


FILE * fp;

// Hann window
void hann(int size);
double *hann_values;

double complex *i_fft_a, *i_time_a, *o_fft_a, *o_time_a;
double complex *i_fft_b, *i_time_b, *o_fft_b, *o_time_b;

fftw_plan i_forward_a, o_inverse_a;
fftw_plan i_forward_b, o_inverse_b;

jack_port_t *input_port, *output_port;
int buffers_size;

jack_client_t *client;


double *freqs;
double sample_rate;

int t_0;
int index_buff_0;

jack_default_audio_sample_t **in, **out,mag_a,mag_b;

jack_port_t **input;
jack_port_t **output;


float mic_distance;

int jack_callback (jack_nframes_t nframes, void *arg)
{
	int i,j,k = 0;
    double ccv[nframes];
	for(i = 0; i < NUM_CH; ++i)
	{
		in[i] = jack_port_get_buffer(input[i],nframes);
		out[i] = jack_port_get_buffer(output[i],nframes);
	}

	for(i = 0; i < nframes; ++i)
	{
		out[0][i] = in[0][i] ;//* hann_values[i];
		out[1][i] = in[1][i] ;//* hann_values[i];  
	}

/*
	for(i = 0; i < nframes; ++i)
    {
    	in[0][i] = exp(-(i)/10.0);	
    }

    for(i = 0; i < nframes; ++i)
    {
    	in[1][i] = 0.0;	
    }

	for(i = 150; i < nframes; ++i)
    {
    	in[1][i] = in[0][i-150];	
    }
*/

	float a_mean = 0;
	float b_mean = 0;


	for(i = 0; i < nframes; ++i)
	{
		in[0][i] = in[0][i] ;//* hann_values[i];
		in[1][i] = in[1][i] ;//* hann_values[i];  
	}

	for(i = 0; i < nframes; ++i)
	{
		a_mean += in[0][i];
		b_mean += in[1][i];
	}

	a_mean /= nframes;
	b_mean /= nframes;


	for(j = 0; j < nframes; ++j)
	{
		i_time_a[j] = in[0][j] - a_mean;
		i_time_b[j] = in[1][j] - b_mean;
	}

	fftw_execute(i_forward_a);
	fftw_execute(i_forward_b);

	for(i = 0; i < nframes; ++i)
	{
		o_fft_a[i] = i_fft_a[i] * conj(i_fft_b[i]) ;	
	}

	fftw_execute(o_inverse_a);

	for(i = 0; i < nframes; ++i)
	{
		mag_a += pow(in[0][i],2); 
		mag_b += pow(in[1][i],2); 
	}

	mag_a = sqrt(mag_a);
	mag_b = sqrt(mag_b);


	for(i = 0; i < nframes; ++i)
	{
		ccv[i] = creal(o_time_a[i]) /( mag_a*mag_b);
	}
	
	
	/*
	printf("--------\n");
	for(i = 0; i < nframes; ++i)
		printf("%f, ",ccv[i]);  
	printf("--------\n");
*/

	float max = 0;
	int index;
	for(i = 0; i < nframes; ++i)
	{
		if( ccv[i] > max )
		{
		 max = ccv[i]; 
		 index = i;
		}
	}

	int desface=0;
	float angle;

	if(index < nframes/2 )
		desface = index;
	else
		desface = index -nframes;

	//printf("max: %f index: %d   desface: %d\n",max,index,desface);
	//exit(0);

	angle =((desface/sample_rate)*(343.2))/mic_distance;

	//printf("angle : %lf \n",  angle );
	//printf("%f \n",  (acos(angle)*180)/M_PI );
	printf("%f\n",  angle );
	fprintf(fp, "%f \n", angle);
	fflush(fp);

	return 0;
}


/**
 * JACK calls this shutdown_callback if the server ever shuts down or
 * decides to disconnect the client.
 */
void jack_shutdown (void *arg){
	exit (1);
}


int main (int argc, char *argv[]) {
	const char *client_name = "desp_freq";
	jack_options_t options = JackNoStartServer;
	jack_status_t status;
	int i,j,k = 0;
	char name_aux[50];



	fp = popen("python simulator_node.py", "w");
    if (fp == NULL) {
        printf("popen error\n");
        exit(1);
    }



	mic_distance = atof(argv[1]);

	//JACK PORTS
	input = (jack_port_t**)malloc(NUM_CH*sizeof(jack_port_t*));
	output = (jack_port_t**)malloc(NUM_CH*sizeof(jack_port_t*));

	// IN OUT CALLBACK
	in  = (jack_default_audio_sample_t** )malloc(NUM_CH * sizeof(jack_default_audio_sample_t*));
	out = (jack_default_audio_sample_t** )malloc(NUM_CH *sizeof(jack_default_audio_sample_t*));

	for(i = 0; i < NUM_CH; ++i)
	{
		in[i] = (jack_default_audio_sample_t*)malloc(NUM_CH * sizeof(jack_default_audio_sample_t));
		out[i] = (jack_default_audio_sample_t*)malloc(NUM_CH * sizeof(jack_default_audio_sample_t));	
	}

	/* open a client connection to the JACK server */
	client = jack_client_open (client_name, options, &status);
	if (client == NULL){
		/* if connection failed, say why */
		printf ("jack_client_open() failed, status = 0x%2.0x\n", status);
		if (status & JackServerFailed) {
			printf ("Unable to connect to JACK server.\n");
		}
		exit (1);
	}

	/* if connection was successful, check if the name we proposed is not in use */
	if (status & JackNameNotUnique){
		client_name = jack_get_client_name(client);
		printf ("Warning: other agent with our name is running, `%s' has been assigned to us.\n", client_name);
	}
	
	/* tell the JACK server to call 'jack_callback()' whenever there is work to be done. */
	jack_set_process_callback (client, jack_callback, 0);
	
	
	/* tell the JACK server to call 'jack_shutdown()' if it ever shuts down,
	   either entirely, or if it just decides to stop calling us. */
	jack_on_shutdown (client, jack_shutdown, 0);
	
	/* display the current sample rate. */
	//printf ("Sample rate: %d\n", jack_get_sample_rate (client));
	//printf ("Window size: %d\n", jack_get_buffer_size (client));
	sample_rate = (double)jack_get_sample_rate(client);
	int nframes = jack_get_buffer_size (client);
	
	hann(nframes);
	

	freqs = (double *) malloc(sizeof(double)*nframes);
	for(int i = 0; i <= nframes/2; ++i )
	{
		freqs[i] = i*(sample_rate/nframes);
		if( i > 0 && i < nframes/2 )
			freqs[nframes -i] =  -1*freqs[i];
	}
	

	


	//preparing FFTW3 buffers
	//preparing FFTW3 buffers
	i_fft_a  = (double complex *) fftw_malloc(sizeof(double complex) * nframes );
	i_time_a = (double complex *) fftw_malloc(sizeof(double complex) * nframes );
	o_fft_a  = (double complex *) fftw_malloc(sizeof(double complex) * nframes );
	o_time_a = (double complex *) fftw_malloc(sizeof(double complex) * nframes );
	
	i_forward_a = fftw_plan_dft_1d(nframes , i_time_a, i_fft_a , FFTW_FORWARD, FFTW_MEASURE);
	o_inverse_a = fftw_plan_dft_1d(nframes , o_fft_a , o_time_a, FFTW_BACKWARD, FFTW_MEASURE);

	i_fft_b  = (double complex *) fftw_malloc(sizeof(double complex) * nframes );
	i_time_b = (double complex *) fftw_malloc(sizeof(double complex) * nframes );
	o_fft_b  = (double complex *) fftw_malloc(sizeof(double complex) * nframes );
	o_time_b = (double complex *) fftw_malloc(sizeof(double complex) * nframes );
	
	i_forward_b = fftw_plan_dft_1d(nframes, i_time_b, i_fft_b , FFTW_FORWARD, FFTW_MEASURE);
	o_inverse_b = fftw_plan_dft_1d(nframes, o_fft_b , o_time_b, FFTW_BACKWARD, FFTW_MEASURE);
	
	/* create the agent input port */
	for(i = 0; i< NUM_CH; ++i)
	{
		sprintf(name_aux,"input%d",i);
		input[i]= jack_port_register(client,name_aux,JACK_DEFAULT_AUDIO_TYPE,JackPortIsInput,0);
		sprintf(name_aux,"output%d",i);
		output[i]= jack_port_register(client,name_aux,JACK_DEFAULT_AUDIO_TYPE,JackPortIsOutput,0);
	}

	for(i = 0; i< NUM_CH; ++i)
		if (input[i] == NULL || output[i] == NULL   ) 
		{
			printf("Could not create agent ports. Have we reached the maximum amount of JACK agent ports?\n");
			exit (1);
		}

	/* Tell the JACK server that we are ready to roll.
	   Our jack_callback() callback will start running now. */

	if (jack_activate (client)) {
		printf ("Cannot activate client.");
		exit (1);
	}
	

	//printf ("Agent activated.\n");
	
	/* Connect the ports.  You can't do this before the client is
	 * activated, because we can't make connections to clients
	 * that aren't running.  Note the confusing (but necessary)
	 * orientation of the driver backend ports: playback ports are
	 * "input" to the backend, and capture ports are "output" from
	 * it.
	 */
	//printf ("Connecting ports... ");
	 
	/* Assign our input port to a server output port*/
	// Find possible output server port names
	const char **serverports_names;


        serverports_names = jack_get_ports (client, NULL, NULL, JackPortIsOutput);
	//serverports_names = jack_get_ports (client, NULL, NULL, JackPortIsPhysical|JackPortIsOutput);
	if (serverports_names == NULL) {
		printf("No available physical capture (server output) ports.\n");
		exit (1);
	}
	// Connect the first available to our input port
	
	for(i=0;i<=3;i++)
		printf("Aqui  %s \n",serverports_names[i]);
	for(i = 1; i <= NUM_CH; ++i)
	{
		if (jack_connect (client, serverports_names[i], jack_port_name (input[i-1]))) {
			printf("Cannot connect input port.\n");
			exit (1);
		}
	}
	// free serverports_names variable for reuse in next part of the code
	free (serverports_names);
	
	
	/* Assign our output port to a server input port*/
	// Find possible input server port names
	serverports_names = jack_get_ports (client, NULL, NULL, JackPortIsPhysical|JackPortIsInput);
	if (serverports_names == NULL) {
		printf("No available physical playback (server input) ports.\n");
		exit (1);
	}
	// Connect the first available to our output port
	for(i = 0; i < NUM_CH; ++i)
	{
		if (jack_connect (client, jack_port_name (output[i]), serverports_names[i])) {
			printf("Cannot connect input port.\n");
			exit (1);
		}
	}
	// free serverports_names variable, we're not going to use it again
	free (serverports_names);
	
	
	//printf ("done.\n");
	/* keep running until stopped by the user */
	sleep (-1);
	
	
	/* this is never reached but if the program
	   had some other way to exit besides being killed,
	   they would be important to call.
	*/
	jack_client_close (client);
	exit (0);
}


void hann(int size )
{
	int i;
	
	hann_values = (double *)malloc(sizeof(double)*size);
	for ( i = 0; i < size; i++) 
	    hann_values[i] = 0.5 * (1 - cos( 2 * M_PI * i/ (size - 1.0) ));
}





