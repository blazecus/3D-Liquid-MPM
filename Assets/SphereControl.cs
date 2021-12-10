using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class SphereControl : MonoBehaviour
{
    // Start is called before the first frame update
    void Start()
    {
        
    }

    // Update is called once per frame
    void Update()
    {
        Rigidbody rb = GetComponent<Rigidbody>();
        if (Input.GetKey(KeyCode.A))
            rb.AddForce(Vector3.left*20);
        else if (Input.GetKey(KeyCode.D))
            rb.AddForce(Vector3.right*20);
        else if (Input.GetKey(KeyCode.W))
            rb.AddForce(Vector3.forward * 20);
        else if (Input.GetKey(KeyCode.S))
            rb.AddForce(Vector3.back * 20);
        else
            rb.velocity = new Vector3(0, 0, 0);
    }
}
