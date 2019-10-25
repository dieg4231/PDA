#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <math.h>

#include <jack/jack.h>

// Include FFTW header
#include <Eigen/Eigen>
#include <Eigen/Eigenvalues> 
#include <complex> //needs to be included before fftw3.h for compatibility
#include <fftw3.h>

#include <iostream>


#define WINDOW_SIZE 1024
#define WINDOWS_PER_BUFF_0 6
#define WINDOWS_PER_BUFF_AB 4
#define NUM_CH  2
#define N_SAMPLES_X_MUSIC 2
#define SOUND_SPEED 343

//using Eigen::MatrixXd;




FILE * fp;

// Hann window
void hann(int size);
double *hann_values;

std::complex<double> *i_fft_a, *i_time_a, *o_fft_a, *o_time_a;
std::complex<double> *i_fft_b, *i_time_b, *o_fft_b, *o_time_b;

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



Eigen::MatrixXcd x(NUM_CH,N_SAMPLES_X_MUSIC*512);
Eigen::MatrixXcd x_ct(N_SAMPLES_X_MUSIC,NUM_CH);
Eigen::MatrixXcd x_aux(NUM_CH,NUM_CH);
Eigen::MatrixXcd st_vec(NUM_CH,180);
std::complex<double> music_spectrum[180];
std::complex<double> s1[1000],s2[1000];

double t[1000];

int cta_n_samples_x_music=0;
int max;
double anlges[180];
double freq =  1500;
std::complex<double> imag(0.0,1.0);
int cta = 0;
int jack_callback (jack_nframes_t nframes, void *arg)
{
	int i,j,k = 0;
	for(i = 0; i < NUM_CH; ++i)
	{
		in[i] = (jack_default_audio_sample_t *)jack_port_get_buffer(input[i],nframes);
		out[i] = (jack_default_audio_sample_t *)jack_port_get_buffer(output[i],nframes);
	}


	for(j = 0; j < nframes; ++j)
	{
		out[0][j] = in[0][j];
		out[1][j] = in[1][j];
	}


	for(j = 0; j < nframes; ++j)
	{
		i_time_a[j] = in[0][j];
		i_time_b[j] = in[1][j];
	}


	fftw_execute(i_forward_a);
	fftw_execute(i_forward_b);
	
	if( cta_n_samples_x_music < N_SAMPLES_X_MUSIC)
	{
		for (i = 0; i < 512; ++i)
		{
		x(0,cta)= (i_fft_a[i]);//1500 hz
		x(1,cta)= (i_fft_b[i]);//1500 hz
		++cta;
		}
		++cta_n_samples_x_music;
	}
	else
	{
		cta = 0;

		freq =  freqs[31];
		std::cout << "___________________" << std::endl;
		std::cout << "xxxxxxxxxxx" << std::endl;
		std::cout << x << std::endl;
		std::cout << "TTTTTTTTTTT" << std::endl;
		x_ct = x.conjugate().transpose();
		std::cout << x_ct << std::endl;
		std::cout << "MULLLLLLLLLLLLLLLLLL" << std::endl;
		std::cout << x*x_ct << std::endl;
		x_aux = (x*x_ct)/N_SAMPLES_X_MUSIC;


      /*
		x_aux(0,0)= std::complex<double>(1.00000 , 0.00000);
		x_aux(0,1)= std::complex<double>(0.32426 , 0.94597);
		
		x_aux(1,0)= std::complex<double>(0.32426 ,- 0.94597);
		x_aux(1,1)= std::complex<double>(1.00000 , 0.00000);


*/
		/*
		x_aux(0,0)= std::complex<double>(2.0 , 0.0000);
		x_aux(0,1)= std::complex<double>( -0.58451 , 1.36328 );
		x_aux(0,2)= std::complex<double>(-0.13800 ,- 0.14499);
		
		x_aux(1,0)= std::complex<double>(-0.58451 ,- 1.36328);
		x_aux(1,1)= std::complex<double>(2.00000 , 0.00000);
		x_aux(1,2)= std::complex<double>(-0.58451 , 1.36328);

		x_aux(2,0)= std::complex<double>(-0.13800 , 0.14499);
		x_aux(2,1)= std::complex<double>(-0.58451 ,- 1.36328);
		x_aux(2,2)= std::complex<double>(2.00000 , 0.00000);
		*/
		std::cout << "XX:\n" << x_aux << std::endl ;


		Eigen::ComplexEigenSolver<Eigen::MatrixXcd>  es(x_aux);
		std::cout << "The eigenvalues of A are:\n" << es.eigenvalues().real() << std::endl << std::endl;
		std::cout << "The eigenvectors of A are (one vector per column):\n" << es.eigenvectors().real() << std::endl << std::endl;


		max =0;// ((double)es.eigenvalues()[0].real() > (double)es.eigenvalues()[1].real() )? 1:0;

		std::cout << "El ruidoso:  " << max << " \n" << std::endl << std::endl;


		for( i = 0; i <180; ++i)
		{
			st_vec(0,i) = 1;
			st_vec(1,i) = exp( -imag * (double)2.0 * M_PI * freq * (double)(mic_distance / SOUND_SPEED) * sin((i-90)*M_PI/180) );
		//	st_vec(2,i) = exp( -imag * (double)2.0 * M_PI * freq * (double)(2.0*mic_distance / SOUND_SPEED) * sin((i-90)*M_PI/180) );	     	
		} 
		
		for( i = 0; i <180; ++i)
		{
			music_spectrum[i]= (std::complex<double>)( st_vec.col(i).transpose()*st_vec.col(i)  )/(std::complex<double>)(st_vec.col(i).transpose()*es.eigenvectors().col(max)*es.eigenvectors().col(max).transpose()*st_vec.col(i) );
			std::cout << abs(music_spectrum[i])   << std::endl;
		}





		cta_n_samples_x_music=0;
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
	const char *client_name = "MUMUMUSIC";
	jack_options_t options = JackNoStartServer;
	jack_status_t status;
	int i,j,k = 0;
	char name_aux[50];


	for( i = 0; i <180; ++i)
		anlges[i] = i;  


	/*fp = popen("python simulator_node.py", "w");
    if (fp == NULL) {
        printf("popen error\n");
        exit(1);
    }*/



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
	i_fft_a  = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * nframes );
	i_time_a = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * nframes );
	o_fft_a  = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * nframes );
	o_time_a = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * nframes );
	
	i_forward_a = fftw_plan_dft_1d(nframes , reinterpret_cast<fftw_complex*>(i_time_a), reinterpret_cast<fftw_complex*>(i_fft_a ), FFTW_FORWARD, FFTW_MEASURE);
	o_inverse_a = fftw_plan_dft_1d(nframes , reinterpret_cast<fftw_complex*>(o_fft_a ), reinterpret_cast<fftw_complex*>(o_time_a), FFTW_BACKWARD, FFTW_MEASURE);

	i_fft_b  = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * nframes );
	i_time_b = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * nframes );
	o_fft_b  = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * nframes );
	o_time_b = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * nframes );
	
	i_forward_b = fftw_plan_dft_1d(nframes, reinterpret_cast<fftw_complex*>(i_time_b), reinterpret_cast<fftw_complex*>(i_fft_b ), FFTW_FORWARD, FFTW_MEASURE);
	o_inverse_b = fftw_plan_dft_1d(nframes, reinterpret_cast<fftw_complex*>(o_fft_b ), reinterpret_cast<fftw_complex*>(o_time_b), FFTW_BACKWARD, FFTW_MEASURE);
	
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

	for(i = 0; i < NUM_CH; ++i)
	{
		if (jack_connect (client, serverports_names[i], jack_port_name (input[i]))) {
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





