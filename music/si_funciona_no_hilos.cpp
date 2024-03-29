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
#define NUM_CH  2
#define N_SAMPLES_X_MUSIC 8
#define SOUND_SPEED 343
#define RANGE 180
#define N_FRECS 4

FILE * fp;

std::complex<double> *i_fft_a, *i_time_a, *o_fft_a, *o_time_a;
std::complex<double> *i_fft_b, *i_time_b, *o_fft_b, *o_time_b;
std::complex<double> *i_fft_c, *i_time_c, *o_fft_c, *o_time_c;

fftw_plan i_forward_a, o_inverse_a;
fftw_plan i_forward_b, o_inverse_b;
fftw_plan i_forward_c, o_inverse_c;

jack_port_t *input_port, *output_port;
jack_client_t *client;

double *freqs;
double sample_rate;

jack_default_audio_sample_t **in, **out;

jack_port_t **input;
jack_port_t **output;


float mic_distance;



Eigen::MatrixXcd *x;//(NUM_CH,N_SAMPLES_X_MUSIC*512) [6];
Eigen::MatrixXcd x_ct(N_SAMPLES_X_MUSIC,NUM_CH);
Eigen::MatrixXcd x_aux(NUM_CH,NUM_CH);
Eigen::MatrixXcd st_vec(NUM_CH,RANGE);
std::complex<double> music_spectrum[N_FRECS][RANGE];



int cta_n_samples_x_music=0;
int max;
double anlges[RANGE];
double freq ;
std::complex<double> imag(0.0,1.0);
int index_frex[] = {4,6,7,8,51,61};

float avg=0;
int cta_callado=0;
double max_energy =0;
std::complex<double> norm;
int max_energy_index = 0;
bool once = true;

int jack_callback (jack_nframes_t nframes, void *arg)
{
	int i,j,k = 0;

	avg=0;

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
	//	i_time_c[j] = in[2][j];	
	}

	fftw_execute(i_forward_a);
	fftw_execute(i_forward_b);
	//fftw_execute(i_forward_c);

	/* Boton rojo */
	for( i = 1; i <= nframes / 2; ++i ) avg += i_fft_a[i].real();
	avg /= nframes / 2;
	avg > 0.001 ? cta_callado=0  : cta_callado ++ ;
	if( cta_callado > 5)
	{
		fprintf(fp,"*\n");
		fflush(fp);
		return 0;
	}
	/*--------------------------------------------*

	/*Las frecuencias más reprecentativas*/
	if(once)
	{

		once = false;
		max_energy =0;

		for( i = 1; i <= nframes/20; ++i )
			if ( max_energy < fabs(i_fft_a[i].real()) )
			{
				index_frex[0] = i;
				max_energy=fabs(i_fft_a[i].real());
			}
		for (j = 1; j < N_FRECS; ++j)
		{
			max_energy = 0;
			for(int i = 0; i <= nframes/20; ++i )
			{
				if ( max_energy < fabs(i_fft_a[i].real()) &&  fabs(i_fft_a[i].real()) <  fabs(i_fft_a[index_frex[j-1]].real()) &&  fabs(i_fft_a[i].real()) !=  fabs(i_fft_a[index_frex[j-1]].real()) )
				{
					index_frex[j] = i;
					max_energy=abs(i_fft_a[i].real());
				}
			}
		}

		std::cout << "Frecuencias seleccionadas -----"  << std::endl;
		for (j = 0; j < N_FRECS; ++j)
			std::cout << "FF " << index_frex[j] << " " <<  freqs[index_frex[j]] << " " << fabs(i_fft_a[index_frex[j]].real()) << std::endl;
	}
	/*---------------------------------------------------------------*/

	if( cta_n_samples_x_music < N_SAMPLES_X_MUSIC)
	{
		for (i = 0; i < N_FRECS; ++i)
		{
			
			x[i](0,cta_n_samples_x_music)= (i_fft_a[index_frex[i]]);
			x[i](1,cta_n_samples_x_music)= (i_fft_b[index_frex[i]]);
			//x[i](2,cta)= (i_fft_c[index_frex[i]]);//1500 hz
		}
		++cta_n_samples_x_music;
	}
	else
	{
		once = true;
		std::cout << "___________________" << std::endl;
		for(j=0;j<N_FRECS;j++)
		{

			freq =  freqs[index_frex[j]];
			x_ct = x[j].transpose().conjugate();
			x_aux = (x[j]*x_ct)/N_SAMPLES_X_MUSIC;

			Eigen::ComplexEigenSolver<Eigen::MatrixXcd>  es(x_aux);
			
			/*Encdontrando  el indice del eigenvalor mas chico*/
			max =0;
			for(i =1; i<NUM_CH;i++)
				if( es.eigenvalues()[i].real() < es.eigenvalues()[max].real()  )
					max = i;
			
			for( i = 0; i < RANGE ; ++i)
			{
				st_vec(0,i) = 1;
				st_vec(1,i) = exp( -imag * (double)2.0 * M_PI * freq * (double)(mic_distance / SOUND_SPEED) *  -sin( ((float)i-90)*M_PI/180.0) );
				//st_vec(2,i) = exp( -imag * (double)2.0 * M_PI * freq * (double)(mic_distance / SOUND_SPEED) * -cos( ((float)(150-i) )*M_PI/180.0) );	     	
				//std::cout <<"vv "<< (float)i*M_PI/180.0 << std::endl;
			}

			for( i = 0; i < RANGE ; ++i)
				music_spectrum[j][i]= (std::complex<double>)( st_vec.col(i).transpose()*st_vec.col(i)  )/(std::complex<double>)(st_vec.col(i).transpose()*es.eigenvectors().col(max)*es.eigenvectors().col(max).transpose()*st_vec.col(i) );
			
			/*Normalizando*/
			/*
			max=0;
			for( i = 1; i < RANGE ; ++i)
				if( music_spectrum[j][max].real() < music_spectrum[j][i].real() )
					max = i;

			norm = music_spectrum[j][max];

			for( i = 0; i < RANGE ; ++i)
				music_spectrum[j][i] /= norm ;
			*/
			
		}

		for( j = 1; j < N_FRECS ; ++j)
			for( i = 0; i < RANGE ; ++i)
				music_spectrum[0][i] *= music_spectrum[j][i];

		for( i = 0; i < RANGE ; ++i)
		{
			fprintf(fp,"%.3f\n", abs(music_spectrum[0][i]));
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

	for( i = 0; i < RANGE ; ++i)
		st_vec(0,i) = 1;
		

	for( i = 0; i <180; ++i)
		anlges[i] = i;  

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
	for(i = 0; i < NUM_CH; ++i)
	{
		/*
		if (jack_connect (client, jack_port_name (output[i]), serverports_names[i])) {
			printf("Cannot connect input port.\n");
			exit (1);
		}
		*/
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

