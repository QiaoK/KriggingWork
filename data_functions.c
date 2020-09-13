/*
* Copyright (C) 2016, Northwestern University.
* This file contains data I/O interface for the rest of program.
* See also datafunctions.h
*/

#include "clusterfunctions.h"
#include "datafunctions.h"
#include "krigfunctions.h"

//Definition of constants
#define STATIONS "//stations.txt"
#define STATION_RECORD_LENGTH 20
#define STATION_META_DATE_LENGTH 81
#define DATA_LENGTH 36
#define HEADER_LENGTH 24
#define X_COORDINATE_LENGTH 8
#define Y_COORDINATE_LENGTH 7
#define ATTRIBUTE_LENGTH 5
#define DATE_LENGTH 10
#define LEVEL_LENGTH 4
#define BUFFER_SIZE 4096
/*
 * This function is called by read_IGRA. It reads a single station.
 * data: The data array to be stored.
 * time_elapse: Number of time stamps
 * cx: Spatial coordinate x.
 * cy: Spatial coordinate y.
 * local_copy: The local file stream to be stored to. The value is NULL is no local copy is required.

*/
void read_IGRA_station(Object** data,char* path,DWORD time_elapse,char* cx,char* cy,FILE* local_copy){
	FILE* stream = fopen(path, "r"); /* should check the result */
	//Initialize strings
	char header[HEADER_LENGTH+2];
	header[HEADER_LENGTH]='\0';
	char record[DATA_LENGTH+1];
	record[DATA_LENGTH]='\0';
	char attribute[ATTRIBUTE_LENGTH+2];
	attribute[ATTRIBUTE_LENGTH]='\0';
	char date[DATE_LENGTH+1];
	char comma=',',end='\n';
	char c_records[LEVEL_LENGTH+1];
	char buffer[BUFFER_SIZE];
	c_records[LEVEL_LENGTH]='\0';
	DWORD i,j,counter;
	DWORD n_records;
	//printf("check\n");
	DTYPE x=atof(cx);
	DTYPE y=atof(cy);
	//Fill in data for each time stamp at this station.
	for(i=0;i<time_elapse;i++){
		fgets(header, HEADER_LENGTH+2,stream);
		//printf("line=%d,last=%c\n",strlen(header),header[HEADER_LENGTH]);
		memcpy(date,header+6,sizeof(char)*DATE_LENGTH);
		//2 dimensional coordinates
		data[i]=Calloc(1,Object);
		data[i]->spatial_coordinates=Calloc(2,SPATIAL_TYPE);
		data[i]->spatial_coordinates[0]=x;
		data[i]->spatial_coordinates[1]=y;
		data[i]->time=atof(date);
		data[i]->neighbors=-1;
		data[i]->normalized_value=0;
		memcpy(c_records,header+20,sizeof(char)*LEVEL_LENGTH);
		n_records=atoi(c_records);
		//printf("c_records=%s\n",c_records);
		counter=0;
		data[i]->attribute=0;
		//Sum up all records
		for(j=0;j<n_records;j++){
			fgets(record, DATA_LENGTH+2, stream);
			//printf("%ld\n",j);
			//printf("temperature type=%c\n",record[20]);
			if(record[20]=='B'&&record[0]=='2'&&record[1]=='1'){
				memcpy(attribute,record+15,sizeof(char)*ATTRIBUTE_LENGTH);
				//printf("temperature=%s\n",temperature);
				if(atof(attribute)!=-9999){
					data[i]->attribute+=atof(attribute);
					counter++;
				}
			}
		}
		//Take average if data exists.
		if(counter!=0){
			data[i]->attribute/=(counter*10);
			strcpy(buffer,"0000000000");
			sprintf(buffer,"%.*lf",2,data[i]->attribute);
			//printf("%ld\n",counter);
		}else{
			//-9999 means this attribute does not exist
			data[i]->attribute=-9999;
			strcpy(buffer,"-9999");
		}
		//Store to local file stream if needed.
		if(local_copy!=NULL){
			fwrite(cx,X_COORDINATE_LENGTH,sizeof(char),local_copy);
			fwrite(&comma,1,sizeof(char),local_copy);
			fwrite(cy,Y_COORDINATE_LENGTH,sizeof(char),local_copy);
			fwrite(&comma,1,sizeof(char),local_copy);
			fwrite(buffer,ATTRIBUTE_LENGTH+1,sizeof(char),local_copy);
			fwrite(&comma,1,sizeof(char),local_copy);
			fwrite(date,DATE_LENGTH,sizeof(char),local_copy);
			fwrite(&end,1,sizeof(char),local_copy);
		}
	}
}

void write_clusters(Clusters* clusters,char* filename,BOOLEAN header){
	printf("start writing clusters to local file system\n");
	DWORD i;
	Objects* objects;
	char buffer[BUFFER_SIZE];
	strcpy(buffer,filename);
	DWORD len=strlen(buffer);
	for(i=0;i<clusters->size;i++){
		objects=cluster_to_objects(clusters->clusters[i]);
		sprintf(buffer+len,"%lld.txt",i);
		write_spatial_temporal_data(objects,buffer,header);
		Free(objects->objects);
		Free(objects);
	}
}

void write_variogram_result(char* filename, BOOLEAN header, Samples* samples, DTYPE* C, VARIOGRAM_TYPE variogram_type){
	FILE* local_copy = fopen( filename , "w" );
	DWORD i;
	DTYPE parameters[2];
	DTYPE gamma_hat;
	if(header){
		fprintf(local_copy,"%s,%s,%s,%s,%s\n","x","phi","N","y","predict");
	}
	for(i=0;i<samples->size;i++){
		parameters[0]=samples->x[i];
		parameters[1]=samples->phi[i];
		gamma_hat=compute_variogram_by_parameters(parameters,C,variogram_type);
		fprintf(local_copy,"%lf,%lf,%lld,%lf,%lf\n",samples->x[i],samples->phi[i],samples->N[i],samples->y[i],gamma_hat);
	}
	fclose(local_copy);
}

void write_spatial_temporal_data(Objects* data,char* filename,BOOLEAN header){
	FILE* local_copy = fopen( filename , "w" );
	DWORD i;
	if(header){
		fprintf(local_copy,"%s,%s,%s,%s\n","x","y","time","attribute");
	}
	for(i=0;i<data->size;i++){
		fprintf(local_copy,"%lf,%lf,%lf,%lf\n",data->objects[i]->spatial_coordinates[0],data->objects[i]->spatial_coordinates[1],data->objects[i]->time,data->objects[i]->attribute);
	}
	fclose(local_copy);
}

Objects* read_IGRA(char* dir_name,DWORD time_elapse,DWORD station_size,BOOLEAN local_copy_flag){
	//Intialize strings with constant size.

	DWORD size=strlen(dir_name);
	char stations[size+strlen(STATIONS)+1];
	//printf("%d\n",strlen(STATIONS));
	strcpy(stations,dir_name);
	strcpy(stations+size,STATIONS);
	//printf("%s\n", stations);
	FILE* file = fopen(stations, "r"); /* open meta file for IGRA dataset */
	char line[STATION_META_DATE_LENGTH+3];
	line[STATION_META_DATE_LENGTH+2]='\0';
	char station_number[128];
	station_number[5]='\0';
	strcpy(station_number+5,".y2d");
	/*
	Example format
	"US_IGRA"+"//"+"72201.y2d//72201.y2d"+"\0",
	"72201.y2d//72201.y2d" has size STATION_LENGTH
	*/
	char* station_record=Calloc(size+STATION_RECORD_LENGTH+3,char);
	char* x=Calloc(X_COORDINATE_LENGTH+1,char);
	char* y=Calloc(Y_COORDINATE_LENGTH+1,char);
	x[X_COORDINATE_LENGTH]='\0';
	y[Y_COORDINATE_LENGTH]='\0';
	DWORD i=0;
	FILE *local_copy=NULL;
	if(local_copy_flag){
		local_copy = fopen( "local_copy.csv" , "w" );
	}
	Object** data=Calloc(station_size*time_elapse,Object*);
	Object** current=data; /*Pointer to the tail of array data.*/
	//Iterate through all stations from meta data file
	while (i<station_size&&fgets(line, sizeof(line), file)) {
		//printf("length=%d,line=%s",strlen(line), line);
		memcpy(station_number,line+4,sizeof(char)*5);
		//printf("%s\n",station_number);
		memcpy(station_record,dir_name,sizeof(char)*size);
		memcpy(station_record+size,"//",sizeof(char)*2);
		memcpy(station_record+size+2,station_number,sizeof(char)*9);
		memcpy(station_record+size+11,"//",sizeof(char)*2);
		strcpy(station_record+size+13,station_number);
		//printf("%s\n",station_number);
		//printf("%s\n",station_record);
		memcpy(x,line+54,sizeof(char)*X_COORDINATE_LENGTH);
		memcpy(y,line+47,sizeof(char)*Y_COORDINATE_LENGTH);
		//printf("x=%s,y=%s\n",x,y);
		//Read data o=frin a statuib,
		if(!access(station_record, F_OK)){
			//printf("%s\n",station_record);
			read_IGRA_station(current,station_record,time_elapse,x,y,local_copy);
			current+=time_elapse;
			i++;
		}
		//printf("longitude=%s,latitude=%s\n",x,y);
	}
	if(i==0){
		/*No data has been read, return NULL*/
		return NULL;
	}
	if(local_copy_flag){
		fclose(local_copy);
	}
	data=(Object**)realloc(data,sizeof(Object*)*i*time_elapse);
	Objects* result=Calloc(1,Objects);
	result->objects=data;
	result->size=i*time_elapse;
	fclose(file);
	Free(station_record);
	Free(x);
	Free(y);
	return result;
}
/*
 * Count number of lines in a string file.
 * filename: The file to be processed.
 * Return: Number of lines in that file.
*/
DWORD count_lines(char* filename){
	DWORD lines=0;
	FILE* fd=fopen(filename,"r");
	char c;
	while((c=fgetc(fd))!=EOF){
		if(c=='\n'){
			lines++;
		}
	}
	return lines;
}
/*
 * Convert unix time format to temporal data type.
 * time: Unix time in string format
 * Return: Seconds since 1970/01/01 0:0
*/
TEMPORAL_TYPE parse_time(char* time){
	struct tm tm;
	strptime(time, "%Y-%m-%d  %H:%M:%S", &tm);
	time_t t = mktime(&tm);
	return (DTYPE) t;
}
/*
 * Convert a spatio-temporal point in string format to Object pointer format.
 * A spatio-temporal object should contain x,y coordinates, a time stamp and a single non-spatial attribute.
 * line: String that contains data for the spatio-temporal point.
 * Return: spatio-temporal point in object pointer format.
*/
Object* read_spatial_temporal_object(char* line){
	Object* object=Calloc(1,Object);
	object->spatial_coordinates=Calloc(2,SPATIAL_TYPE);
	//parse for x coordinate
	char* cpt=line;
	DWORD size=0;
	while(*cpt!=','){
		cpt++;
		size++;
	}
	cpt++;
	char x[size+1];
	x[size]='\0';
	memcpy(x,line,sizeof(char)*size);
	object->spatial_coordinates[0]=atof(x);
	//parse for y coordinate
	line=cpt;
	size=0;
	while(*cpt!=','){
		cpt++;
		size++;
	}
	cpt++;
	char y[size+1];
	y[size]='\0';
	memcpy(y,line,sizeof(char)*size);
	object->spatial_coordinates[1]=atof(y);
	//parse for time value
	line=cpt;
	size=0;
	while(*cpt!=','){
		cpt++;
		size++;
	}
	cpt++;
	char time[size+1];
	time[size]='\0';
	memcpy(time,line,sizeof(char)*size);
	//printf("time to be processed=%s\n",time);
	object->time=parse_time(time);
	//parse for attribute
	size=0;
	line=cpt;
	while(*cpt!='\0'){
		cpt++;
		size++;
	}
	cpt++;
	char t[size+1];
	t[size]='\0';
	memcpy(t,line,sizeof(char)*size);
	object->attribute=atof(t);
	return object;
}
/*
 * Convert a spatial point in string format to Object pointer format.
 * A spatio-temporal object should contain x,y coordinatesand a single non-spatial attribute.
 * line: String that contains data for the spatial point.
 * Return: Spatial point in object pointer format.
*/
Object* read_object(char* line){
	Object* object=Calloc(1,Object);
	object->spatial_coordinates=Calloc(2,SPATIAL_TYPE);
	char* cpt=line;
	//locate coma
	DWORD size=0;
	while(*cpt!=','){
		cpt++;
		size++;
	}
	cpt++;
	//Read x coordinate
	char x[size+1];
	x[size]='\0';
	memcpy(x,line,sizeof(char)*size);
	line=cpt;
	object->spatial_coordinates[0]=atof(x);
	size=0;
	while(*cpt!=','){
		cpt++;
		size++;
	}
	cpt++;
	//Read y coordinate
	char y[size+1];
	y[size]='\0';
	memcpy(y,line,sizeof(char)*size);
	object->spatial_coordinates[1]=atof(y);
	size=0;
	line=cpt;
	while(*cpt!='\0'){
		cpt++;
		size++;
	}
	cpt++;
	//Read attribute
	char t[size+1];
	t[size]='\0';
	memcpy(t,line,sizeof(char)*size);
	object->attribute=atof(t);
	return object;
}

Objects* read_csv(char* filename,WORD type,BOOLEAN header){
	/*Double stream reading algorithm*/
	FILE* stream1=fopen(filename,"r");
	FILE* stream2=fopen(filename,"r");
	DWORD lines=count_lines(filename);
	Object** objects=Calloc(lines,Object*);
	DWORD i=0;
	DWORD length=0;
	char* line;
	char c;
	/*Escape headers if there are any*/
	if(header){
		fgetc(stream2);
		while((c=fgetc(stream1))!='\n'){
			fgetc(stream2);
		}
	}
	while((c=fgetc(stream1))!=EOF){
		length++;
		if(c=='\n'&&length>=1){
			//Use the second stream to obtain a line.
			line=Calloc(length+1,char);
			fgets(line,length,stream2);
			//printf("%s\n",line);
			fgetc(stream2);
			line[length]='\0';
			//Parse the line to object pointer
			switch(type){
				case SPATIAL_DATA:
					objects[i]=read_object(line);
					break;
				case SPATIAL_TEMPORAL_DATA:
					objects[i]=read_spatial_temporal_object(line);
			}
			//Initialize the object
			objects[i]->neighbors=-1;
			objects[i]->normalized_value=0;
			i++;
			length=0;
			Free(line);
		}
	}
	//Clean up and return
	fclose(stream1);
	fclose(stream2);
	Objects* result=Calloc(1,Objects);
	result->size=i;
	result->objects=objects;
	return result;
}

/*
int main(void){
	Objects* data=read_IGRA("US_IGRA",60,221,TRUE);
	//print_data(data);
}

void print_data(Objects* data){
	DWORD i;
	for(i=0;i<data->size;i++){
		printf("x=%lf,y=%lf,attribute=%lf,time=%lf\n",data->objects[i]->spatial_coordinates[0],data->objects[i]->spatial_coordinates[1],data->objects[i]->attribute,data->objects[i]->time);
	}
}

void main(){
	Objects* data=read_csv("pollution_data.csv",SPATIAL_TEMPORAL_DATA,TRUE);
	print_data(data);
	//printf("%lf\n",parse_time("2011-12-12 22:59:48"));
	//printf("%lf\n",parse_time("2011-12-12 22:52:48"));
}
*/
