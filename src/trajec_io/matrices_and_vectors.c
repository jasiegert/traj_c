void matrix33_cofactors(float mat[3][3], float adj[3][3])
{
    adj[0][0] = mat[1][1] * mat[2][2] - mat[2][1] * mat[1][2];
    adj[0][1] = mat[1][0] * mat[2][2] - mat[2][0] * mat[1][2];
    adj[0][2] = mat[1][0] * mat[2][1] - mat[2][0] * mat[1][1];
    adj[1][0] = mat[0][1] * mat[2][2] - mat[2][1] * mat[0][2];
    adj[1][1] = mat[0][0] * mat[2][2] - mat[2][0] * mat[0][2];
    adj[1][2] = mat[0][0] * mat[2][1] - mat[2][0] * mat[0][1];
    adj[2][0] = mat[0][1] * mat[1][2] - mat[1][1] * mat[0][2];
    adj[2][1] = mat[0][0] * mat[1][2] - mat[1][0] * mat[0][2];
    adj[2][2] = mat[0][0] * mat[1][1] - mat[1][0] * mat[0][1];
    return;
}

void matrix33_transpose(float mat[3][3], float transpose[3][3])
{
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            transpose[i][j] = mat[j][i];
        }
    }
    return;
}

void matrix33_multiplication(float a[3][3], float b[3][3], float c[3][3])
{
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            c[i][j] = 0;
            for (int k = 0; k < 3; k++)
            {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    return;
}

void matrix33_inverse(float mat[3][3], float inv[3][3])
{
    float cofactors[3][3];
    matrix33_cofactors(mat, cofactors);
    matrix33_transpose(cofactors, inv);
    float determinant = 0;
    for (int i = 0; i < 3; i++)
    {
        determinant += mat[0][i] * cofactors[0][i];
    }
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            inv[i][j] /= determinant;
        }
    }
    return;
}

void matrix33_vector3_multiplication(float mat[3][3], float vec[3], float out[3])
{
    for (int i = 0; i < 3; i++)
    {
        out[i] = 0;
        for (int j = 0; j < 3; j++)
        {
            out[i] += mat[i][j] * vec[j];
        }
    }
    return;
}

