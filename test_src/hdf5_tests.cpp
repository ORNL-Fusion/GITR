#include <iostream>
#include <vector>
#include <string>
#include <cassert>
#include "catch2/catch_all.hpp"
#include "hdf5.h"
#include "interpolator.h"

TEST_CASE( "hdf5" )
{
  SECTION( "t0" )
  {  
    int rank = 3;

    hid_t file;

    std::string file_name = "captain.h5";
    std::string dataset_name = "/data";

    file = H5Fcreate( file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );

    hsize_t current_dims[ 3 ] = { 10, 11, 13};

    hsize_t max_dims[ 3 ] = { 10, 11, 13 };

    double data_buffer_0[ 10 ][ 11 ][ 13 ];

    double count = 0;
    for( int i = 0; i < 10; i++ )
    {
      for( int ii = 0; ii < 11; ii++ )
      {
        for( int iii = 0; iii < 13; iii++ )
        {
          data_buffer_0[ i ][ ii ][ iii ] = count++;
        }
      }
    }

    hid_t space_id = H5Screate( H5S_NULL );

    H5Sset_extent_simple( space_id, rank, current_dims, max_dims );

    /* create modifiable datatype - you can change the endianness if you wanted here */
    hid_t data_type = H5Tcopy( H5T_NATIVE_DOUBLE );

    double fill_value = 0;

    hid_t plist = H5Pcreate( H5P_DATASET_CREATE );

    herr_t status = H5Pset_fill_value( plist, data_type, &fill_value );

    hid_t data_set = H5Dcreate( file, dataset_name.c_str(), 
                                data_type, space_id, H5P_DEFAULT, plist, H5P_DEFAULT );

    status = H5Dwrite( data_set, data_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_buffer_0 );

    /* close literally everything */
    H5Dclose(data_set);
    
    H5Sclose( space_id );

    H5Tclose( data_type );

    H5Pclose( plist );
    /* end close */

    /* open the same file again...  */
    unsigned int mode = H5F_ACC_RDWR;
    hid_t file_id = H5Fopen( file_name.c_str(), mode, H5P_DEFAULT );
    REQUIRE( file_id != H5I_INVALID_HID );

    /* Captain! Put it in a c_str()! */
    hid_t dataset_id = H5Dopen( file_id, dataset_name.c_str(), H5P_DEFAULT );
    REQUIRE( dataset_id != H5I_INVALID_HID );

    /* get the dataspace: Captain! This does not even appear to be necessary...  */
    hid_t dataspace_id = H5Dget_space( dataset_id );
    REQUIRE( dataspace_id != H5I_INVALID_HID );

    /* this is the rank */
    int const n_dataspace_dims = H5Sget_simple_extent_ndims( dataspace_id );
    REQUIRE( n_dataspace_dims > -1 );

    /* get the dimensions of the dataset, name must be enforced */
    std::vector< hsize_t > dims( n_dataspace_dims );
    std::vector< hsize_t > dummy( n_dataspace_dims );
    int tmp = H5Sget_simple_extent_dims( dataspace_id, dims.data(), dummy.data() );
    REQUIRE( tmp == n_dataspace_dims );

    int size = H5Sget_simple_extent_npoints( dataspace_id );
    REQUIRE( size > -1 );

    /* create a buffer for the data to be read into: it will obviously be of size npoints  */
    /* can this be one of your multidimensional arrays?? Why not?? Just do it! */

    hid_t datatype_id = H5Dget_type( dataset_id );
    REQUIRE( datatype_id != H5I_INVALID_HID );

    /* according to storage format, this data maps to: leading dimension first, zyx */
    double data_buffer_1[ 10 ][ 11 ][ 13 ];
       
    /* read the data back into a different buffer */

    status = H5Dread( dataset_id, datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_buffer_1 );
    REQUIRE( status > -1 );

    /* create a flat tensor out of data_buffer_1 and the dimensions */
    tensor< double > t( (double *)(data_buffer_1), dims.data(), n_dataspace_dims );

    int index = 0;
    /* print out the contents of the file */
    for( int i = 0; i < 10; i++ )
    {
      for( int ii = 0; ii < 11; ii++ )
      {
        for( int iii = 0; iii < 13; iii++ )
        {
          long long unsigned int point[ 3 ] = { i, ii, iii };
          /* Captain! Make the buffer multidimensional */
          REQUIRE( t.get( point ) == data_buffer_1[ i ][ ii ][ iii ] );
        }
      }
    }

    /* Captain! Next we have to pull in the interpolator tests - copy those below */
  }
}





























