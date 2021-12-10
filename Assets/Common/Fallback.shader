Shader "Custom/Fallback" {
    SubShader{
      Tags { "RenderType" = "Opaque" }
      CGPROGRAM
      #pragma surface surf Lambert addshadow
      struct Input {
          float4 color : COLOR;
      };
      void surf(Input IN, inout SurfaceOutput o) {
          o.Albedo = 1;
      }
      ENDCG
    }
        Fallback "Diffuse"
}
