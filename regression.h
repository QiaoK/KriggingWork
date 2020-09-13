
#ifndef KRIG_REGRESSION_H
#define KRIG_REGRESSION_H
extern Matrix* construct_st_model(Objects* data,DWORD time_lag,DWORD neighbor_number);
extern Objects* getNearestNeighbors(Object* target,Objects* data,DWORD number);
extern Objects* locate_map(Objects* data,DTYPE time);
extern Objects* get_previous_time_stamps(Object* target,Objects* data,DWORD number);
extern DTYPE predict_attribute(Objects* data,Object* target,DWORD time_lag,DWORD neighbor_number,Matrix* coefficients);

#endif
