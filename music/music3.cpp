#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <math.h>

#include <jack/jack.h>

// Include FFTW header
#include <thread> 
#include <Eigen/Eigen>
#include <Eigen/Eigenvalues> 
#include <complex> //needs to be included before fftw3.h for compatibility
#include <fftw3.h>

#include <iostream>

#define WINDOW_SIZE 1024
#define WINDOWS_PER_BUFF_0 6
#define WINDOWS_PER_BUFF_AB 4
#define NUM_CH  3
#define N_SAMPLES_X_MUSIC 6
#define SOUND_SPEED 343.0
#define RANGE 360
#define N_FRECS 3

//using Eigen::MatrixXd;

FILE * fp;

// Hann window
void hann(int size);
double *hann_values;

std::complex<double> *i_fft_a, *i_time_a, *o_fft_a, *o_time_a;
std::complex<double> *i_fft_b, *i_time_b, *o_fft_b, *o_time_b;
std::complex<double> *i_fft_c, *i_time_c, *o_fft_c, *o_time_c;

fftw_plan i_forward_a, o_inverse_a;
fftw_plan i_forward_b, o_inverse_b;
fftw_plan i_forward_c, o_inverse_c;

jack_port_t *input_port, *output_port;
int buffers_size;

jack_client_t *client;


double *freqs;
double sample_rate;


jack_default_audio_sample_t **in, **out,mag_a,mag_b;

jack_port_t **input;
jack_port_t **output;


float mic_distance;

Eigen::MatrixXcd *x;//(NUM_CH,N_SAMPLES_X_MUSIC) [6];
std::complex<double>  **music_spectrum;

int cta_n_samples_x_music=0;

std::complex<double> imag(0.0,1.0);
int index_frex[] = {31,101,321,131,51,61};


void music(double freq,Eigen::MatrixXcd x,std::complex<double> * music_spectrumt)
{
	Eigen::MatrixXcd x_ct(N_SAMPLES_X_MUSIC,NUM_CH);
	Eigen::MatrixXcd x_aux(NUM_CH,NUM_CH);
	Eigen::MatrixXcd st_vec(NUM_CH,RANGE);
	int min =0;
	int i;

	x_ct = x.conjugate().transpose();

	x_aux = (x*x_ct)/N_SAMPLES_X_MUSIC;

	Eigen::ComplexEigenSolver<Eigen::MatrixXcd>  es(x_aux);
	
	for(i = 1; i < 3; i++)
		if( es.eigenvalues()[i].real() < es.eigenvalues()[min].real()  )
			min = i;
	
	for( i = 0; i < RANGE ; ++i)
	{
		st_vec(0,i) = 1;
		st_vec(1,i) = exp( -imag * (double)2.0 * M_PI * freq * (double)(mic_distance / SOUND_SPEED) *  sin( ((float)i-180)*M_PI/180.0) );
		st_vec(2,i) = exp( -imag * (double)2.0 * M_PI * freq * (double)(mic_distance / SOUND_SPEED) * -cos( ((float)i-180)*M_PI/180.0) );	     	
	}
		
	for( i = 0; i < RANGE ; ++i)
		music_spectrumt[i] = abs((std::complex<double>)( st_vec.col(i).transpose()*st_vec.col(i)  )/(std::complex<double>)(st_vec.col(i).transpose()*es.eigenvectors().col(min)*es.eigenvectors().col(min).transpose()*st_vec.col(i) ));
}


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

	if( cta_n_samples_x_music < N_SAMPLES_X_MUSIC)
	{
		for(j = 0; j < nframes; ++j)
		{
			i_time_a[j] = in[0][j];
			i_time_b[j] = in[1][j];
			i_time_c[j] = in[2][j];
		}

		fftw_execute(i_forward_a);
		fftw_execute(i_forward_b);
		fftw_execute(i_forward_c);

		for (i = 0; i < N_FRECS; ++i)
		{
			x[i](0,cta_n_samples_x_music)= (i_fft_a[index_frex[i]]);//1500 hz
			x[i](1,cta_n_samples_x_music)= (i_fft_b[index_frex[i]]);//1500 hz
			x[i](2,cta_n_samples_x_music)= (i_fft_c[index_frex[i]]);//1500 hz		
		}
		++cta_n_samples_x_music;
	}
	else
	{
		std::thread th0 (music,index_frex[0],x[0],music_spectrum[0]);
		std::thread th1 (music,index_frex[1],x[1],music_spectrum[1]); 
		std::thread th2 (music,index_frex[2],x[2],music_spectrum[2]); 
		//std::thread th3 (music,index_frex[3],x[3],music_spectrum[3]); 
		//std::thread th4 (music,index_frex[4],x[4],music_spectrum[4]); 
		//std::thread th5 (music,index_frex[5],x[5],music_spectrum[5]); 

		th0.join(); 
		th1.join(); 
		th2.join(); 
		//th3.join(); 
		//th4.join(); 
		//th5.join(); 
		
		std::cout << "-----"  << std::endl;

		for( j = 1; j < N_FRECS ; ++j)
			for( i = 0; i < RANGE ; ++i)
				music_spectrum[0][i] += music_spectrum[j][i];

		//std::cout << std::fixed;
    	//std::cout << std::setprecision(2);
		for( i = 0; i < RANGE ; ++i)
		//	std::cout << music_spectrum[0][i].real()  << std::endl;
		{
			fprintf(fp,"%.3f\n", music_spectrum[0][i].real());
			fflush(fp);
		}
		fprintf(fp,"-----\n");
		fflush(fp);
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


	fp = popen("python simulator_node.py", "w");
    if (fp == NULL) {
        printf("popen error\n");
        exit(1);
    }

	mic_distance = atof(argv[1]);


	x = (Eigen::MatrixXcd *)malloc(N_FRECS*sizeof(Eigen::MatrixXcd));

	for(i=0;i<N_FRECS;i++)
	{
		x[i]= Eigen::MatrixXcd(NUM_CH,N_SAMPLES_X_MUSIC) ;
	}

	music_spectrum = (std::complex<double>  **)malloc(N_FRECS*sizeof(std::complex<double>  *));
	
	for(i = 0; i < N_FRECS; i++)
		music_spectrum[i] = (std::complex<double> *) malloc(RANGE*sizeof(std::complex<double> ));

	

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

	i_fft_c  = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * nframes );
	i_time_c = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * nframes );
	o_fft_c  = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * nframes );
	o_time_c = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * nframes );
	
	i_forward_c = fftw_plan_dft_1d(nframes, reinterpret_cast<fftw_complex*>(i_time_c), reinterpret_cast<fftw_complex*>(i_fft_c ), FFTW_FORWARD, FFTW_MEASURE);
	o_inverse_c = fftw_plan_dft_1d(nframes, reinterpret_cast<fftw_complex*>(o_fft_c ), reinterpret_cast<fftw_complex*>(o_time_c), FFTW_BACKWARD, FFTW_MEASURE);
	
	/* create the agent input port */
	for(i = 0; i< NUM_CH; ++i)
	{
		sprintf(name_aux,"input%d",i+1);
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
		//if (jack_connect (client, serverports_names[i], jack_port_name (input[i]))) {
		//	printf("Cannot connect input port.\n");
		//	exit (1);
		//}
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
	for(i = 0; i < 2; ++i)
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





