/**
 * Desface en frecuencia 

      t0  t1  t2  t3  t4  t5
buff0 |---|---|---|---|---|---|

      t5  t0  t1  t2  t3  t4  
buff0 |---|---|---|---|---|---|

      t4  t5  t0  t1  t2  t3    
buff0 |---|---|---|---|---|---|

      t3  t4  t5  t0  t1  t2      
buff0 |---|---|---|---|---|---|

+++++++++++++++++++++++++++++++++++++++++++++++

      t0  t1  t2  t3  t4  t5
buff0 |---|---|---|---|---|---|


  t0  t1  t2  t3  
A |---|---|---|---|

		  t2  t3  t4  t5
        B |---|---|---|---|
            |   |
            |   |
            |   |
            |---|
            
  |--2560---|                
(1024+1024+512)
 
 2560 / 48000Hz = 0.053 s



 **/

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <jack/jack.h>

// Include FFTW header
#include <complex.h> //needs to be included before fftw3.h for compatibility
#include <fftw3.h>

#define WINDOW_SIZE 1024
#define WINDOWS_PER_BUFF_0 6
#define WINDOWS_PER_BUFF_AB 4
#define NUM_CH  2
// Hann window
void hann(int size);
double *hann_values;

double complex *i_fft_a, *i_time_a, *o_fft_a, *o_time_a,**buffer_0;
double complex *i_fft_b, *i_time_b, *o_fft_b, *o_time_b;
fftw_plan i_forward_a, o_inverse_a;
fftw_plan i_forward_b, o_inverse_b;

jack_port_t *input_port, *output_port;
int buffers_size;

jack_client_t *client;


double *freqs;
double sample_rate;

double desface;

int t_0;
int index_buff_0;

jack_default_audio_sample_t **in, **out;

jack_port_t **input;
jack_port_t **output;


int jack_callback (jack_nframes_t nframes, void *arg)
{
	int i,j,k = 0;

	//printf("NFRAMES %d\n",nframes );
	for(i = 0; i < NUM_CH; ++i)
	{
		in[i] = jack_port_get_buffer(input[i],nframes);
		out[i] = jack_port_get_buffer(output[i],nframes);
	}
	// Corrimiento de los apuntadores del buffer 0
	for(i = 0; i < WINDOWS_PER_BUFF_0 - 1; ++i )
		buffer_0[i] = buffer_0[i + 1];

	buffer_0[5] = buffer_0[0];

	//nuevos datos
	for (int i = 0; i < WINDOW_SIZE; ++i)
		buffer_0[5][i] = in[0][i];

	// LLenamos buffer a y b con la señal hanneada
	for(i = 0, k = 0; i < WINDOWS_PER_BUFF_AB; ++i)
		for(j = 0; j < WINDOW_SIZE; ++j, ++k)
			{
				i_time_a[k] = buffer_0[i][j] * hann_values[k];
				i_time_b[k] = buffer_0[i+2][j] * hann_values[k];
			}

	fftw_execute(i_forward_a);
	fftw_execute(i_forward_b);

	//Desplazamiento
	for(i = 0; i < WINDOW_SIZE * WINDOWS_PER_BUFF_AB; ++i)
	{
		o_fft_a[i] = i_fft_a[i] * cexp(-I * 2 * M_PI * freqs[i] * desface);
		o_fft_b[i] = i_fft_b[i] * cexp(-I * 2 * M_PI * freqs[i] * desface);	
	}

	fftw_execute(o_inverse_a);
	fftw_execute(o_inverse_b);

	// sumamos señales

    for(i = WINDOW_SIZE / 2, j =  WINDOW_SIZE * WINDOWS_PER_BUFF_AB - (3 * WINDOW_SIZE / 2), k = 0 ; k < WINDOW_SIZE; ++i, ++j, ++k)
    {
    	out[0][k] = (creal(o_time_b[i]) + creal(o_time_a[j]) )/ (WINDOW_SIZE * WINDOWS_PER_BUFF_AB);	

    }

    //Señal sin dezplazamiento
    
    for(i = WINDOW_SIZE / 2, k = 0 ; i < WINDOW_SIZE; ++i, ++k)
    {
    	out[1][k] = creal(buffer_0[2][i]);	
    }

    for(i = 0  ; i < WINDOW_SIZE / 2; ++i, ++k)
    {
    	out[1][k] = creal(buffer_0[3][i]);	
    }

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



	//JACK PORTS
	input = (jack_port_t**)malloc(NUM_CH*sizeof(jack_port_t*));
	output = (jack_port_t**)malloc(NUM_CH*sizeof(jack_port_t*));

	//BUFF DESPLAZAMIENTOS
	buffer_0 = ( double complex **) malloc( sizeof(double complex *) * WINDOWS_PER_BUFF_0 );
	for(i = 0; i < WINDOWS_PER_BUFF_0 ; ++i)
		buffer_0[i] = (double complex*)malloc(sizeof(double complex) * WINDOW_SIZE);	

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
	printf ("Sample rate: %d\n", jack_get_sample_rate (client));
	printf ("Window size: %d\n", jack_get_buffer_size (client));
	sample_rate = (double)jack_get_sample_rate(client);
	int nframes = jack_get_buffer_size (client);

	if(nframes != WINDOW_SIZE)
	{
		printf("Sorry nframes and WINDOW_SIZE  must be the same size. :( \n");
	}

	if( argc > 1)
	{
		desface = atof(argv[1]);
		if( fabs(desface) > (nframes/jack_get_sample_rate (client)) )
		{
			printf("Parameters ok\n");
		}
		else
		{
			printf(" %lf  OUT OF RANGE !!! \n Range values ( -%f  , %f ) \n",desface,(float)(nframes/jack_get_sample_rate (client)),(float)(nframes/jack_get_sample_rate (client)));
			return 0;
		}
	}
	else
	{
		printf(" :( Please provide a number -.- \n");
		return 0;
	}

	
	hann(nframes*4);

	
	/*
	for(int i = 0 ; i < 100; i++)
		printf("%f\n",hann_values[i]);
	return 0;
	*/

	//
	freqs = (double *) malloc(sizeof(double)*nframes* WINDOWS_PER_BUFF_AB);
	for(int i = 0; i <= (nframes* WINDOWS_PER_BUFF_AB)/2; ++i )
	{
		freqs[i] = i*(sample_rate/(nframes* WINDOWS_PER_BUFF_AB));
		if( i > 0 && i < (nframes* WINDOWS_PER_BUFF_AB)/2 )
			freqs[(nframes* WINDOWS_PER_BUFF_AB) -i] =  -1*freqs[i];
	}
	
	//preparing FFTW3 buffers
	i_fft_a  = (double complex *) fftw_malloc(sizeof(double complex) * nframes * WINDOWS_PER_BUFF_AB);
	i_time_a = (double complex *) fftw_malloc(sizeof(double complex) * nframes * WINDOWS_PER_BUFF_AB);
	o_fft_a  = (double complex *) fftw_malloc(sizeof(double complex) * nframes * WINDOWS_PER_BUFF_AB);
	o_time_a = (double complex *) fftw_malloc(sizeof(double complex) * nframes * WINDOWS_PER_BUFF_AB);
	
	i_forward_a = fftw_plan_dft_1d(nframes * WINDOWS_PER_BUFF_AB, i_time_a, i_fft_a , FFTW_FORWARD, FFTW_MEASURE);
	o_inverse_a = fftw_plan_dft_1d(nframes * WINDOWS_PER_BUFF_AB, o_fft_a , o_time_a, FFTW_BACKWARD, FFTW_MEASURE);

	i_fft_b  = (double complex *) fftw_malloc(sizeof(double complex) * nframes * WINDOWS_PER_BUFF_AB);
	i_time_b = (double complex *) fftw_malloc(sizeof(double complex) * nframes * WINDOWS_PER_BUFF_AB);
	o_fft_b  = (double complex *) fftw_malloc(sizeof(double complex) * nframes * WINDOWS_PER_BUFF_AB);
	o_time_b = (double complex *) fftw_malloc(sizeof(double complex) * nframes * WINDOWS_PER_BUFF_AB);
	
	i_forward_b = fftw_plan_dft_1d(nframes * WINDOWS_PER_BUFF_AB, i_time_b, i_fft_b , FFTW_FORWARD, FFTW_MEASURE);
	o_inverse_b = fftw_plan_dft_1d(nframes * WINDOWS_PER_BUFF_AB, o_fft_b , o_time_b, FFTW_BACKWARD, FFTW_MEASURE);
	
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
	

	printf ("Agent activated.\n");
	
	/* Connect the ports.  You can't do this before the client is
	 * activated, because we can't make connections to clients
	 * that aren't running.  Note the confusing (but necessary)
	 * orientation of the driver backend ports: playback ports are
	 * "input" to the backend, and capture ports are "output" from
	 * it.
	 */
	printf ("Connecting ports... ");
	 
	/* Assign our input port to a server output port*/
	// Find possible output server port names
	const char **serverports_names;
	serverports_names = jack_get_ports (client, NULL, NULL, JackPortIsPhysical|JackPortIsOutput);
	if (serverports_names == NULL) {
		printf("No available physical capture (server output) ports.\n");
		exit (1);
	}
	// Connect the first available to our input port
	for(i = 0; i < NUM_CH; ++i)
	{
		if (jack_connect (client, serverports_names[0], jack_port_name (input[i]))) {
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
	
	
	printf ("done.\n");
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





