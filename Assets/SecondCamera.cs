using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class SecondCamera : MonoBehaviour
{
    [SerializeField] private Shader depthShader;
    [SerializeField] private string replacementTag;
    private RenderTexture renderTex;
    private void OnEnable()
    {
        var cam = GetComponent<Camera>();
        if (cam == null)
            cam = GetComponent<Camera>();

        renderTex = new RenderTexture(cam.pixelWidth, cam.pixelHeight, 24);
        renderTex.antiAliasing = Mathf.Max(1, QualitySettings.antiAliasing);

        cam.targetTexture = renderTex;

        if (depthShader != null)
            cam.SetReplacementShader(depthShader, replacementTag);

        Shader.SetGlobalTexture("_CustomDepthTexture", renderTex);
    }
}
